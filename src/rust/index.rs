use anyhow::{Context, Result};
use needletail::{parse_fastx_file, parse_fastx_stdin};
use rayon::prelude::*;
use rocksdb::{Options, DB};
use serde::{Deserialize, Serialize};

use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::str;
use std::sync::{Arc, Mutex};

use crate::aminoacid::AminoAcidAmbiguity;
use crate::encoding::{get_encoding_fn_from_moltype, get_hash_function_from_moltype};
use crate::kmer_signature::{KmerInfo, SerializableSignature, SignatureKmerMapping};
use crate::protein::Protein;
use crate::uniprot::UniProtEntry;
use sourmash::_hash_murmur;
use sourmash::cmd::ComputeParameters;
use sourmash::collection::Collection;
use sourmash::manifest::Manifest;
use sourmash::signature::Signature;
use sourmash::storage::{FSStorage, InnerStorage};
use sourmash_plugin_branchwater::search_significance::{
    compute_inverse_document_frequency, get_hash_frequencies, Normalization,
};
use sourmash_plugin_branchwater::utils::multicollection::SmallSignature;

/// Statistics for k-mer frequency analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteomeIndexKmerStats {
    pub idf: f64,       // Inverse document frequency
    pub frequency: f64, // Raw frequency
}

// Represents the serializable state of ProteomeIndex
#[derive(Serialize, Deserialize)]
struct ProteomeIndexState {
    // Store signatures instead of Collection
    signatures: Vec<SerializableSignature>,
    combined_minhash: KmerMinHash,
    protein_info: HashMap<String, SignatureKmerMapping>,
    moltype: String,
    ksize: u32,
    scaled: u32,
}

pub struct ProteomeIndex {
    // RocksDB instance for persistent storage
    db: DB,

    // Collection of all signatures for searching
    collection: Arc<Mutex<Collection>>,

    // Combined minhash of all proteins for statistics
    combined_minhash: Arc<Mutex<KmerMinHash>>,

    // Map of signature md5 -> protein k-mer information
    protein_info: Arc<Mutex<HashMap<String, Protein>>>,

    // Amino acid ambiguity handler
    aa_ambiguity: Arc<AminoAcidAmbiguity>,

    // Protein encoding function
    encoding_fn: fn(u8) -> u8,

    // Add moltype field for serialization
    moltype: String,

    // Add ksize field for serialization
    ksize: u32,

    // Add scaled field for serialization
    scaled: u32,
}

impl ProteomeIndex {
    pub fn new<P: AsRef<Path>>(path: P, ksize: u32, scaled: u32, moltype: &str) -> Result<Self> {
        // Create RocksDB options
        let mut opts = Options::default();
        opts.create_if_missing(true);

        // Open the database
        let db = DB::open(&opts, path)?;

        let hash_function = match get_hash_function_from_moltype(moltype) {
            Ok(value) => value,
            Err(value) => return value,
        };

        let encoding_fn = match get_encoding_fn_from_moltype(moltype) {
            Ok(value) => value,
            Err(value) => return value,
        };

        // Create the minhash sketch
        let minhash = KmerMinHash::new(
            scaled,
            ksize,
            hash_function,
            42,   // seed
            true, // track_abundance
            0,    // num (use scaled instead)
        );

        // Create an empty collection with storage
        let manifest = Manifest::default();
        let storage = InnerStorage::new(
            FSStorage::builder()
                .fullpath("".into())
                .subdir("".into())
                .build(),
        );
        let collection = Collection::new(manifest, storage);

        Ok(Self {
            db,
            collection: Arc::new(Mutex::new(collection)),
            combined_minhash: Arc::new(Mutex::new(minhash)),
            protein_info: Arc::new(Mutex::new(HashMap::new())),
            aa_ambiguity: Arc::new(AminoAcidAmbiguity::new()),
            encoding_fn,
            moltype: moltype.to_string(),
            ksize,
            scaled,
        })
    }

    pub fn add_uniprot_entry(&mut self, uniprot_entry: UniProtEntry) -> Result<()> {
        // Get parameters from combined_minhash
        let mut minhash_guard = self.combined_minhash.lock().unwrap();
        let ksize = minhash_guard.ksize() as u32;
        let scaled = minhash_guard.scaled();

        // Create compute parameters for this protein's signature
        let params = ComputeParameters::builder()
            .ksizes(vec![ksize])
            .scaled(scaled)
            .protein(true)
            .num_hashes(0) // use scaled instead
            .build();

        // Create a signature for this protein
        let mut sig = Signature::from_params(&params);
        sig.set_name(&uniprot_entry.id);

        // Convert to SmallSignature
        let small_sig = SmallSignature {
            location: sig.filename().clone(),
            name: sig.name().clone().unwrap(),
            md5sum: sig.md5sum().to_string(),
            minhash: sig.minhash().unwrap().clone(),
        };
        let md5sum = small_sig.md5sum.clone();
        let mut minhash = small_sig.minhash.clone();
        minhash.add_sequence(uniprot_entry.sequence.as_bytes(), true)?;

        // Add the newly created hashes to the combined minhash
        minhash_guard.add_many_with_abund(&minhash.to_vec_abunds())?;
        drop(minhash_guard);

        // Process k-mers for this protein
        let protein_kmers = self.process_protein_kmers(&uniprot_entry.sequence, &small_sig)?;

        let protein = Protein::new(uniprot_entry, protein_kmers);

        // Update global structures
        {
            let mut collection = self.collection.lock().unwrap();
            let sigs = vec![sig];
            let new_collection = Collection::from_sigs(sigs)?;
            *collection = new_collection;

            let mut proteins = self.protein_info.lock().unwrap();
            proteins.insert(md5sum.clone(), protein.clone());
        }

        // Store in RocksDB
        self.db.put(
            format!("protein:{}", md5sum).as_bytes(),
            bincode::serialize(&protein)?,
        )?;

        Ok(())
    }

    pub fn process_protein_kmers(
        &self,
        sequence: &str,
        signature: &SmallSignature,
    ) -> Result<SignatureKmerMapping> {
        let minhash_guard = signature.minhash.clone();
        let ksize = minhash_guard.ksize() as usize;
        let seed = minhash_guard.seed();
        let hashvals = minhash_guard.to_vec();
        drop(minhash_guard);

        let mut signature_kmers = SignatureKmerMapping::new(SmallSignature {
            location: signature.location.clone(),
            name: signature.name.clone(),
            md5sum: signature.md5sum.clone(),
            minhash: signature.minhash.clone(),
        });

        for i in 0..sequence.len().saturating_sub(ksize - 1) {
            let kmer = &sequence[i..i + ksize];

            // Process the k-mer to get encoded version
            let (encoded_kmer, original_kmer) = self.encode_kmer(kmer);

            // Get the hash from the minhash implementation
            let hashval = _hash_murmur(encoded_kmer.as_bytes(), seed);

            // If this hashval is in the minhash, then save its k-mer positions
            if !hashvals.contains(&hashval) {
                let kmer_info = signature_kmers
                    .kmer_infos
                    .entry(hashval)
                    .or_insert_with(|| KmerInfo {
                        hashval,
                        encoded_kmer: encoded_kmer.clone(),
                        original_kmer_to_position: HashMap::new(),
                    });

                kmer_info
                    .original_kmer_to_position
                    .entry(original_kmer)
                    .or_insert_with(Vec::new)
                    .push(i);
            }
        }

        Ok(signature_kmers)
    }

    /// Process a k-mer to get its encoded version
    pub fn encode_kmer(&self, kmer: &str) -> (String, String) {
        // Create the encoded k-mer by applying the encoding function
        let encoded: Vec<u8> = kmer.chars().map(|c| (self.encoding_fn)(c as u8)).collect();
        let encoded_str = str::from_utf8(&encoded).unwrap();

        (encoded_str.to_owned(), kmer.to_owned())
    }

    /// Get statistics for a hash value
    pub fn get_kmer_stats(&self, hash: u64) -> Result<Option<ProteomeIndexKmerStats>> {
        let key = format!("stats:{}", hash);
        if let Some(value) = self.db.get(key.as_bytes())? {
            Ok(Some(bincode::deserialize(&value)?))
        } else {
            Ok(None)
        }
    }

    /// Process a list of protein FASTA files
    pub fn process_protein_files<P: AsRef<Path> + Sync + Send>(
        &mut self,
        paths: &[P],
    ) -> Result<()> {
        paths.par_iter().try_for_each(|path| -> Result<()> {
            // Create a FASTX reader from the file or stdin
            let mut fastx_reader = if path.as_ref() == Path::new("-") {
                parse_fastx_stdin().context("Failed to parse FASTA/FASTQ data from stdin")?
            } else {
                parse_fastx_file(path).context("Failed to open file for FASTA/FASTQ data")?
            };

            let mut record_count: u64 = 0;

            // Parse records and add sequences to signatures
            while let Some(record_result) = fastx_reader.next() {
                let record = record_result.context("Failed to read a record from input")?;
                let id = String::from_utf8_lossy(record.id()).to_string();
                let seq_bytes = record.seq();
                let sequence = String::from_utf8_lossy(seq_bytes.as_ref()).to_string();

                let protein = UniProtEntry {
                    id,
                    sequence: sequence.clone(),
                    features: Vec::new(),
                    ..Default::default()
                };

                self.process_sequence(seq_bytes.as_ref(), protein)?;
                record_count += 1;
            }
            println!("Processed {} sequence records", record_count);
            Ok(())
        })?;
        Ok(())
    }

    /// Compute and store k-mer statistics
    pub fn compute_statistics(&self) -> Result<()> {
        let minhash = self.combined_minhash.lock().unwrap();
        let collection = self.collection.lock().unwrap();

        // Get all signatures from the collection
        let mut signatures = Vec::new();
        for (_, record) in collection.iter() {
            if let Ok(sig) = collection.sig_from_record(record) {
                let serialized = SerializableSignature::from(sig);
                let small_sig = SmallSignature {
                    location: serialized.location,
                    name: serialized.name,
                    md5sum: serialized.md5sum,
                    minhash: serialized.minhash,
                };
                signatures.push(small_sig);
            }
        }

        // Compute hash frequencies using L1 normalization
        let frequencies = get_hash_frequencies(&minhash, Some(Normalization::L1));

        // Compute inverse document frequency with smoothing
        let idf = compute_inverse_document_frequency(&minhash, &signatures, Some(true));

        // Store statistics for each hash in parallel
        minhash.mins().par_iter().try_for_each(|&hash| {
            let stats = ProteomeIndexKmerStats {
                idf: *idf.get(&hash).unwrap_or(&0.0),
                frequency: *frequencies.get(&hash).unwrap_or(&0.0),
            };

            let key = format!("stats:{}", hash);
            let value = bincode::serialize(&stats)?;
            self.db.put(key.as_bytes(), value)?;
            Ok::<(), anyhow::Error>(())
        })?;

        Ok(())
    }

    /// Process a single protein sequence
    pub fn process_sequence(&self, sequence: &[u8], protein: UniProtEntry) -> Result<()> {
        // Validate sequence - check for invalid amino acids
        for (i, c) in sequence.iter().enumerate() {
            if !self.aa_ambiguity.is_valid_aa(*c as char) {
                return Err(anyhow::anyhow!(
                    "Invalid amino acid '{}' at position {}",
                    *c as char,
                    i + 1
                ));
            }
        }

        // Create a temporary mutable copy of self to add the protein
        let mut opts = Options::default();
        opts.create_if_missing(true);
        let mut temp_self = Self {
            db: DB::open(&opts, self.db.path())?,
            collection: self.collection.clone(),
            combined_minhash: self.combined_minhash.clone(),
            protein_info: self.protein_info.clone(),
            aa_ambiguity: self.aa_ambiguity.clone(),
            encoding_fn: self.encoding_fn,
            moltype: self.moltype.clone(),
            ksize: self.ksize,
            scaled: self.scaled,
        };

        // Add the protein
        temp_self.add_uniprot_entry(protein)?;

        Ok(())
    }

    /// Get all hashes from the combined minhash
    pub fn get_hashes(&self) -> Vec<u64> {
        let minhash = self.combined_minhash.lock().unwrap();
        minhash.mins().to_vec()
    }

    /// Process a UniProt XML file
    pub fn process_uniprot_xml<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut parser = quick_xml::Reader::from_reader(reader);
        let mut buf = Vec::new();

        let mut current_protein = None;
        let mut current_sequence = String::new();

        loop {
            match parser.read_event_into(&mut buf) {
                Ok(quick_xml::events::Event::Start(ref e)) => match e.name().as_ref() {
                    b"entry" => {
                        current_protein = Some(UniProtEntry::default());
                    }
                    b"sequence" => {
                        current_sequence.clear();
                    }
                    _ => {}
                },
                Ok(quick_xml::events::Event::Text(e)) => {
                    if let Some(_protein) = current_protein.as_mut() {
                        if !current_sequence.is_empty() {
                            current_sequence.push_str(&e.unescape()?.into_owned());
                        }
                    }
                }
                Ok(quick_xml::events::Event::End(ref e)) => match e.name().as_ref() {
                    b"entry" => {
                        if let Some(mut protein) = current_protein.take() {
                            protein.sequence = current_sequence.clone();
                            self.process_sequence(current_sequence.as_bytes(), protein)?;
                        }
                        current_sequence.clear();
                    }
                    _ => {}
                },
                Ok(quick_xml::events::Event::Eof) => break,
                Err(e) => return Err(anyhow::anyhow!("Error parsing XML: {}", e)),
                _ => {}
            }
        }
        Ok(())
    }
}
