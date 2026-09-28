use crate::hash_functions::get_hash_function_from_moltype;
use crate::signature::StableSignature;
use crate::types::{KmerSize, MolType};
use crate::SEED;
use serde::{Deserialize, Serialize};
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use std::collections::{HashMap, HashSet};

pub const PROTEIN_TO_MINHASH_RATIO: u32 = 3;

/// ProteinSketches contain Signatures and k-mer position maps
#[derive(Debug, Clone)]
pub struct ProteinSketch {
    name: String,
    signature: StableSignature,
    moltype: MolType,
    protein_ksize: u32,
    scaled: u32,
    // Hashval -> sorted list of positions in the original protein sequence
    kmer_positions: HashMap<u64, Vec<usize>>,
    // Efficient storage data (optional, for performance)
    efficient_data: Option<ProteinSketchStore>,
    // Whether add_protein() should drop low-complexity (homopolymer) k-mers.
    // Defaults to false (legacy behavior); not persisted, since it only affects
    // insertion, not the sketch that results from it.
    remove_low_complexity: bool,
    // Counts from the most recent add_protein() with removal on: k-mer windows
    // examined, and windows removed as low-complexity. Both stay 0 when off.
    // Not persisted; used only for index-time reporting.
    kmer_windows_examined: usize,
    low_complexity_kmers_removed: usize,
}

// Custom serialization for ProteinSketch
impl Serialize for ProteinSketch {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use serde::ser::SerializeStruct;
        let mut state = serializer.serialize_struct("ProteinSketch", 6)?;
        state.serialize_field("name", &self.name)?;
        state.serialize_field("signature", &self.signature)?;
        state.serialize_field("moltype", &self.moltype)?;
        state.serialize_field("protein_ksize", &self.protein_ksize)?;
        state.serialize_field("scaled", &self.scaled)?;
        state.serialize_field("kmer_positions", &self.kmer_positions)?;
        state.end()
    }
}

impl<'de> Deserialize<'de> for ProteinSketch {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use serde::de::{self, MapAccess, Visitor};
        use std::fmt;

        struct ProteinSketchVisitor;

        impl<'de> Visitor<'de> for ProteinSketchVisitor {
            type Value = ProteinSketch;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("struct ProteinSketch")
            }

            fn visit_map<V>(self, mut map: V) -> Result<ProteinSketch, V::Error>
            where
                V: MapAccess<'de>,
            {
                let mut name = None;
                let mut signature = None;
                let mut moltype = None;
                let mut protein_ksize = None;
                let mut kmer_positions = None;
                let mut scaled = None;

                while let Some(key) = map.next_key()? {
                    match key {
                        "name" => {
                            if name.is_some() {
                                return Err(de::Error::duplicate_field("name"));
                            }
                            name = Some(map.next_value()?);
                        }
                        "signature" => {
                            if signature.is_some() {
                                return Err(de::Error::duplicate_field("signature"));
                            }
                            signature = Some(map.next_value()?);
                        }
                        "moltype" => {
                            if moltype.is_some() {
                                return Err(de::Error::duplicate_field("moltype"));
                            }
                            moltype = Some(map.next_value()?);
                        }
                        "protein_ksize" => {
                            if protein_ksize.is_some() {
                                return Err(de::Error::duplicate_field("protein_ksize"));
                            }
                            protein_ksize = Some(map.next_value()?);
                        }
                        "scaled" => {
                            if scaled.is_some() {
                                return Err(de::Error::duplicate_field("scaled"));
                            }
                            scaled = Some(map.next_value()?);
                        }
                        "kmer_positions" => {
                            if kmer_positions.is_some() {
                                return Err(de::Error::duplicate_field("kmer_positions"));
                            }
                            kmer_positions = Some(map.next_value()?);
                        }
                        _ => {
                            map.next_value::<de::IgnoredAny>()?;
                        }
                    }
                }

                let name = name.ok_or_else(|| de::Error::missing_field("name"))?;
                let signature = signature.ok_or_else(|| de::Error::missing_field("signature"))?;
                let moltype = moltype.ok_or_else(|| de::Error::missing_field("moltype"))?;
                let protein_ksize =
                    protein_ksize.ok_or_else(|| de::Error::missing_field("protein_ksize"))?;
                let scaled = scaled.ok_or_else(|| de::Error::missing_field("scaled"))?;
                let kmer_positions =
                    kmer_positions.ok_or_else(|| de::Error::missing_field("kmer_positions"))?;

                Ok(ProteinSketch {
                    name,
                    signature,
                    moltype,
                    protein_ksize,
                    scaled,
                    kmer_positions,
                    efficient_data: None,
                    remove_low_complexity: false,
                    kmer_windows_examined: 0,
                    low_complexity_kmers_removed: 0,
                })
            }
        }

        const FIELDS: &[&str] =
            &["name", "signature", "moltype", "protein_ksize", "scaled", "kmer_positions"];
        deserializer.deserialize_struct("ProteinSketch", FIELDS, ProteinSketchVisitor)
    }
}

/// One residue to the symbol sourmash hashes it as under `moltype`.
///
/// WHY the two cases differ in capitalization: for a table-backed alphabet the sequence
/// is pre-encoded and handed to sourmash as protein, and sourmash uppercases protein
/// input before hashing, so the symbols must be uppercased here to match. For
/// sourmash-encoded alphabets (protein20, dayhoff6, hp_lehninger2) sourmash applies the
/// encoder itself and hashes its lowercase output, so these must NOT be uppercased.
fn hash_encoder(moltype: &str) -> anyhow::Result<impl Fn(u8) -> u8> {
    use crate::alphabets::alphabet_table;
    use crate::hash_functions::get_encoding_fn_from_moltype;

    let residue_classes = alphabet_table(moltype);
    let encoding_fn = get_encoding_fn_from_moltype(moltype)?;
    Ok(move |residue: u8| match residue_classes {
        Some(table) => table
            .get(&residue.to_ascii_uppercase())
            .copied()
            .unwrap_or(residue)
            .to_ascii_uppercase(),
        None => encoding_fn(residue),
    })
}

/// The k-mer windows of one sequence, each expanded into its encoded readings: one for a
/// window with no ambiguous residue, otherwise one per way of resolving it that the
/// alphabet tells apart (see `disambiguate_kmer`, which also drops a window carrying more
/// than `MAX_AMBIGUOUS_RESIDUES_PER_KMER`).
struct Readings<'a> {
    residues: &'a [u8],
    ksize: usize,
    // Checked once for the whole sequence: almost none carry an ambiguous residue (146 of
    // Swiss-Prot 2026_03's 575_748 sequences), so the per-window check is skipped
    // entirely for nearly every sequence.
    has_ambiguous: bool,
}

impl<'a> Readings<'a> {
    fn new(sequence: &'a str, ksize: usize) -> Self {
        use crate::aminoacid::has_ambiguous_residues;

        let residues = sequence.as_bytes();
        Self { residues, ksize, has_ambiguous: has_ambiguous_residues(residues) }
    }

    fn window_count(&self) -> usize {
        (self.residues.len() + 1).saturating_sub(self.ksize)
    }

    /// Call `f` with each window's start position and one encoded reading of it. A
    /// callback rather than an iterator so the unambiguous case, nearly every window, is
    /// encoded into one buffer allocated per sequence instead of copied.
    fn for_each(&self, encode: &impl Fn(u8) -> u8, mut f: impl FnMut(usize, &[u8])) {
        use crate::aminoacid::{disambiguate_kmer, has_ambiguous_residues};

        let mut buffer = Vec::with_capacity(self.ksize);
        for i in 0..self.window_count() {
            let window = &self.residues[i..i + self.ksize];
            if !self.has_ambiguous || !has_ambiguous_residues(window) {
                buffer.clear();
                buffer.extend(window.iter().map(|&residue| encode(residue)));
                f(i, &buffer);
                continue;
            }
            for reading in disambiguate_kmer(window, encode).unwrap_or_default() {
                f(i, &reading);
            }
        }
    }
}

impl ProteinSketch {
    /// # Errors
    ///
    /// Returns an error for an unsupported `moltype`, or for a `protein_ksize`
    /// outside the range `KmerSize` accepts.
    pub fn new(name: &str, protein_ksize: u32, scaled: u32, moltype: &str) -> anyhow::Result<Self> {
        // WHY normalize before anything else: Sourmash-style spellings (`hp`, `dayhoff`)
        // have to pick the same hash function and the same stored name. Reading
        // the hash function from the raw string while storing the normalized one built the
        // minhash for one encoding and pre-encoded for another, which left kmer_positions
        // empty.
        let moltype = MolType::new(moltype).map_err(|e| anyhow::anyhow!(e))?;
        let moltype_str = moltype.get().to_string();
        let hash_function = get_hash_function_from_moltype(&moltype_str)?;
        KmerSize::new(protein_ksize).map_err(|e| anyhow::anyhow!("Invalid k-mer size: {e}"))?;
        let minhash_ksize = protein_ksize * PROTEIN_TO_MINHASH_RATIO;

        let minhash = KmerMinHash::new(
            scaled,
            minhash_ksize,
            hash_function,
            SEED,
            true, // track_abundance
            0,    // num (use scaled instead)
        );

        let signature = StableSignature {
            location: String::new(),
            name: name.to_string(),
            md5sum: String::new(),
            minhash,
            moltype: moltype_str,
            ksize: protein_ksize,
        };

        Ok(Self {
            name: name.to_string(),
            signature,
            moltype,
            protein_ksize,
            scaled,
            kmer_positions: HashMap::new(),
            efficient_data: None,
            remove_low_complexity: false,
            kmer_windows_examined: 0,
            low_complexity_kmers_removed: 0,
        })
    }

    /// Enable or disable dropping low-complexity (homopolymer) k-mers during
    /// `add_protein`. Defaults to `false`.
    pub fn set_remove_low_complexity(&mut self, remove_low_complexity: bool) {
        self.remove_low_complexity = remove_low_complexity;
    }

    /// K-mer windows examined by the most recent `add_protein` with removal on,
    /// and how many of those were removed as low-complexity. Both are 0 when
    /// removal is off, since that path never walks windows itself.
    pub fn low_complexity_counts(&self) -> (usize, usize) {
        (self.kmer_windows_examined, self.low_complexity_kmers_removed)
    }

    /// Create a ProteinSketch from a protein sequence
    pub fn from_protein_sequence(
        name: &str,
        sequence: &str,
        protein_ksize: u32,
        scaled: u32,
        moltype: &str,
    ) -> anyhow::Result<Self> {
        let mut sketch = Self::new(name, protein_ksize, scaled, moltype)?;
        sketch.add_protein(sequence, true)?;
        Ok(sketch)
    }

    /// Create a ProteinSketch from existing signature data
    pub fn from_existing_data(
        name: &str,
        signature: StableSignature,
        moltype: String,
        protein_ksize: u32,
        scaled: u32,
        kmer_positions: HashMap<u64, Vec<usize>>,
    ) -> Self {
        Self {
            name: name.to_string(),
            signature,
            moltype: MolType::new(&moltype).unwrap(),
            protein_ksize,
            scaled,
            kmer_positions,
            efficient_data: None,
            remove_low_complexity: false,
            kmer_windows_examined: 0,
            low_complexity_kmers_removed: 0,
        }
    }

    /// Create a ProteinSketch from efficient storage data
    pub fn from_efficient_data(
        data: ProteinSketchStore,
        moltype: String,
        protein_ksize: u32,
        scaled: u32,
    ) -> anyhow::Result<Self> {
        let hash_function = get_hash_function_from_moltype(&moltype)?;
        let minhash_ksize = protein_ksize * PROTEIN_TO_MINHASH_RATIO;

        let mut minhash = KmerMinHash::new(scaled, minhash_ksize, hash_function, SEED, true, 0);

        if let Some(abunds) = &data.abunds {
            minhash.add_many_with_abund(
                &data.mins.clone().into_iter().zip(abunds.iter().cloned()).collect::<Vec<_>>(),
            )?;
        } else {
            minhash.add_many(&data.mins)?;
        }

        let md5sum = data.mins.iter().fold(0u64, |acc, &min| acc.wrapping_add(min));
        let md5sum = format!("{:x}", md5sum);
        let name = data.name.clone();

        let signature = StableSignature {
            location: String::new(),
            name: name.clone(),
            md5sum,
            minhash,
            moltype: moltype.clone(),
            ksize: protein_ksize,
        };

        Ok(Self {
            name,
            signature,
            moltype: MolType::new(&moltype).unwrap(),
            protein_ksize,
            scaled,
            kmer_positions: data.kmer_positions.clone(),
            efficient_data: Some(data),
            remove_low_complexity: false,
            kmer_windows_examined: 0,
            low_complexity_kmers_removed: 0,
        })
    }

    /// Convert to efficient storage format
    pub fn to_efficient_data(&self, include_raw_sequence: bool) -> ProteinSketchStore {
        let minhash = &self.signature.minhash;
        let mins = minhash.mins().to_vec();
        let abunds = minhash.abunds().map(|abunds| abunds.to_vec());

        let raw_sequence = if include_raw_sequence {
            if let Some(ref data) = self.efficient_data {
                data.raw_sequence.clone()
            } else {
                None
            }
        } else {
            None
        };

        let encoded_sequence = if let Some(ref data) = self.efficient_data {
            data.encoded_sequence.clone()
        } else {
            None
        };

        ProteinSketchStore::new(
            self.signature.name.clone(),
            mins,
            abunds,
            self.kmer_positions.clone(),
            raw_sequence,
            encoded_sequence,
        )
    }

    /// Convert to efficient storage format with pre-allocated sequence capacity
    pub fn to_efficient_data_with_capacity(&self, sequence_capacity: usize) -> ProteinSketchStore {
        let minhash = &self.signature.minhash;
        let mins = minhash.mins().to_vec();
        let abunds = minhash.abunds().map(|abunds| abunds.to_vec());

        ProteinSketchStore::with_sequence_capacity(
            self.signature.name.clone(),
            mins,
            abunds,
            self.kmer_positions.clone(),
            sequence_capacity,
        )
    }

    pub fn set_efficient_data(&mut self, data: ProteinSketchStore) {
        self.efficient_data = Some(data);
    }

    pub fn get_efficient_data(&self) -> Option<&ProteinSketchStore> {
        self.efficient_data.as_ref()
    }

    pub fn has_efficient_data(&self) -> bool {
        self.efficient_data.is_some()
    }

    pub fn get_raw_sequence(&self) -> Option<&str> {
        self.efficient_data.as_ref()?.get_raw_sequence()
    }

    pub fn get_moltype_sequence(&self) -> Option<&str> {
        self.efficient_data.as_ref()?.get_encoded_sequence()
    }

    /// Add a protein sequence, building the minhash and k-mer position map.
    ///
    /// Approach 3: stores `HashMap<u64, Vec<usize>>` (hash → positions) instead of the
    /// old `HashMap<u64, KmerInfo>` (hash → {encoded_kmer, HashMap<orig_kmer, positions>}).
    /// This is ~3.4× faster to build and ~2.5× smaller to serialize, with identical
    /// search speed (O(1) lookup in find_matched_regions).
    pub fn add_protein(&mut self, sequence: &str, store_sequences: bool) -> anyhow::Result<()> {
        use crate::alphabets::alphabet_table;
        use crate::aminoacid::{encode_sequence, has_ambiguous_residues};
        use crate::hash_functions::residue_encoder;

        let moltype_str = self.moltype.to_string();
        let residue_classes = alphabet_table(&moltype_str);
        let encode = hash_encoder(&moltype_str)?;

        // WHY: low-complexity k-mers carry little discriminative signal, so when opted
        // in, a reading that encodes to a run of one class is dropped before insertion.
        // This means hashing one k-mer at a time instead of delegating to sourmash's
        // black-box `add_protein`, which windows and inserts unconditionally.
        // `remove_low_complexity` defaults to false, and the branches below keep every
        // k-mer.
        if self.remove_low_complexity {
            self.add_windows_without_low_complexity(sequence, &encode);
        } else if has_ambiguous_residues(sequence.as_bytes()) {
            // B, J and Z each stand for two residues (B is Asp or Asn, J is Ile or Leu, Z
            // is Glu or Gln). Rather than committing to one, index every window under every
            // reading the alphabet tells apart, so a query carrying either residue matches.
            // sourmash's add_protein windows and hashes internally and cannot do this, so
            // hash window by window here. This branch has to come before the table-backed
            // one, which would otherwise pre-encode the whole sequence in one go and never
            // disambiguate.
            self.add_windows(sequence, &encode);
        } else if let Some(table) = residue_classes {
            // Pre-encode with our custom table so sourmash hashes the reduced symbols via
            // Murmur64Protein (identity). Unknown bytes pass through unchanged.
            let pre_encoded: String = sequence
                .bytes()
                .map(|b| table.get(&b.to_ascii_uppercase()).copied().unwrap_or(b) as char)
                .collect();
            self.signature.minhash.add_protein(pre_encoded.as_bytes())?;
        } else {
            self.signature.minhash.add_protein(sequence.as_bytes())?;
        }

        let md5sum =
            self.signature.minhash.mins().iter().fold(0u64, |acc, &min| acc.wrapping_add(min));
        self.signature.md5sum = format!("{:x}", md5sum);

        self.record_positions(sequence, &encode);

        if store_sequences {
            let efficient_data = self.to_efficient_data_with_capacity(sequence.len());
            let mut efficient_data_with_sequence = efficient_data;
            efficient_data_with_sequence.set_raw_sequence(sequence.to_string());

            let moltype_str = self.moltype.to_string();
            // The full alphabet encodes to itself, so storing an "encoded" copy would just
            // duplicate the raw sequence. MolType normalizes `protein`/`raw` to `protein20`,
            // so this one name covers all three spellings.
            if moltype_str != "protein20" {
                // Display-only, so the table's lowercase symbols are kept as they are, matching
                // built-in hp/dayhoff (which sourmash encodes lowercase); hashing uppercases
                // separately because sourmash uppercases protein input before hashing. An
                // ambiguous residue the alphabet merges is written as its class, one it keeps
                // apart stays as its letter (see encode_sequence).
                let encoded_sequence = String::from_utf8(encode_sequence(
                    sequence.as_bytes(),
                    residue_encoder(&moltype_str),
                ))
                .expect("class symbols and residues are ASCII");
                efficient_data_with_sequence.set_encoded_sequence(encoded_sequence);
            }

            self.set_efficient_data(efficient_data_with_sequence);
        } else {
            let efficient_data = self.to_efficient_data_with_capacity(0);
            self.set_efficient_data(efficient_data);
        }

        Ok(())
    }

    /// Hash every window of the sequence into the minhash, a window carrying an ambiguous
    /// residue under each of its readings.
    fn add_windows(&mut self, sequence: &str, encode: &impl Fn(u8) -> u8) {
        use sourmash::_hash_murmur;

        Readings::new(sequence, self.protein_ksize as usize).for_each(encode, |_, reading| {
            self.signature.minhash.add_hash(_hash_murmur(reading, SEED));
        });
    }

    /// Like `add_windows`, but drops each reading that is low complexity, meaning it
    /// encodes to a run of one class: a raw homopolymer such as `EEEEE` under any
    /// alphabet, and under a reduced alphabet also a window of different residues in one
    /// class, such as `LIVMA` (`hhhhh` under the Lehninger split) or `EEEDD` (`ccccc`
    /// under dayhoff6). Counts the windows walked and the readings dropped for
    /// index-time reporting.
    ///
    /// Readings are expanded here exactly as in `add_windows` and `record_positions`.
    /// Hashing a window with the literal byte B, J or Z instead would put a hash in the
    /// sketch that no canonical query k-mer matches and that the position map, which
    /// looks up each reading's hash, never records.
    fn add_windows_without_low_complexity(&mut self, sequence: &str, encode: &impl Fn(u8) -> u8) {
        use crate::kmer::is_homopolymer_kmer;
        use sourmash::_hash_murmur;

        let readings = Readings::new(sequence, self.protein_ksize as usize);
        self.kmer_windows_examined = readings.window_count();
        self.low_complexity_kmers_removed = 0;
        readings.for_each(encode, |_, reading| {
            if is_homopolymer_kmer(reading) {
                self.low_complexity_kmers_removed += 1;
            } else {
                self.signature.minhash.add_hash(_hash_murmur(reading, SEED));
            }
        });
    }

    /// Record where each minhash min occurs in the sequence, walking the same readings as
    /// the branches that built the minhash so every window hashes to a value they inserted.
    fn record_positions(&mut self, sequence: &str, encode: &impl Fn(u8) -> u8) {
        use sourmash::_hash_murmur;

        let hashvals: HashSet<u64> = self.signature.minhash.mins().iter().copied().collect();
        Readings::new(sequence, self.protein_ksize as usize).for_each(encode, |i, reading| {
            let hashval = _hash_murmur(reading, SEED);
            if hashvals.contains(&hashval) {
                self.kmer_positions.entry(hashval).or_default().push(i);
            }
        });
    }

    pub fn protein_ksize(&self) -> u32 {
        self.protein_ksize
    }

    pub fn moltype(&self) -> &MolType {
        &self.moltype
    }

    pub fn minhash_ksize(&self) -> u32 {
        self.protein_ksize * PROTEIN_TO_MINHASH_RATIO
    }

    pub fn scaled(&self) -> u32 {
        self.scaled
    }

    pub fn into_signature(self) -> StableSignature {
        self.signature
    }

    pub fn signature(&self) -> &StableSignature {
        &self.signature
    }

    /// Get the k-mer positions map: hashval → Vec of positions in the original sequence
    pub fn kmer_positions(&self) -> &HashMap<u64, Vec<usize>> {
        &self.kmer_positions
    }

    /// Get a mutable reference to the k-mer positions map
    pub fn kmer_positions_mut(&mut self) -> &mut HashMap<u64, Vec<usize>> {
        &mut self.kmer_positions
    }

    pub fn is_compatible(&self, other: &ProteinSketch) -> bool {
        self.moltype() == other.moltype()
            && self.scaled() == other.scaled()
            && self.protein_ksize() == other.protein_ksize()
    }

    pub fn mins_as_set(&self) -> HashSet<u64> {
        self.signature().minhash.mins().iter().cloned().collect()
    }

    pub fn intersect(&self, other: &ProteinSketch) -> HashSet<u64> {
        let self_mins = self.mins_as_set();
        let other_mins = other.mins_as_set();
        self_mins.intersection(&other_mins).cloned().collect()
    }
}

/// Efficient storage structure for protein sketches
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(crate = "serde")]
pub struct ProteinSketchStore {
    pub name: String,
    pub mins: Vec<u64>,
    pub abunds: Option<Vec<u64>>,
    /// Hash → positions in the original protein sequence
    pub kmer_positions: HashMap<u64, Vec<usize>>,
    pub raw_sequence: Option<String>,
    pub encoded_sequence: Option<String>,
}

impl ProteinSketchStore {
    pub fn new(
        name: String,
        mins: Vec<u64>,
        abunds: Option<Vec<u64>>,
        kmer_positions: HashMap<u64, Vec<usize>>,
        raw_sequence: Option<String>,
        encoded_sequence: Option<String>,
    ) -> Self {
        Self { name, mins, abunds, kmer_positions, raw_sequence, encoded_sequence }
    }

    pub fn with_sequence_capacity(
        name: String,
        mins: Vec<u64>,
        abunds: Option<Vec<u64>>,
        kmer_positions: HashMap<u64, Vec<usize>>,
        sequence_capacity: usize,
    ) -> Self {
        let raw_sequence = if sequence_capacity > 0 {
            Some(String::with_capacity(sequence_capacity))
        } else {
            None
        };
        Self { name, mins, abunds, kmer_positions, raw_sequence, encoded_sequence: None }
    }

    pub fn set_raw_sequence(&mut self, sequence: String) {
        if self.raw_sequence.is_some() {
            self.raw_sequence = Some(sequence);
        }
    }

    pub fn set_encoded_sequence(&mut self, sequence: String) {
        self.encoded_sequence = Some(sequence);
    }

    pub fn get_raw_sequence(&self) -> Option<&str> {
        self.raw_sequence.as_deref()
    }

    pub fn get_encoded_sequence(&self) -> Option<&str> {
        self.encoded_sequence.as_deref()
    }

    pub fn has_raw_sequence_storage(&self) -> bool {
        self.raw_sequence.is_some()
    }

    pub fn kmer_count(&self) -> usize {
        self.mins.len()
    }

    pub fn estimated_size(&self) -> usize {
        let mut size = self.name.len();
        size += self.mins.len() * 8;
        if let Some(ref abunds) = self.abunds {
            size += abunds.len() * 8;
        }
        for positions in self.kmer_positions.values() {
            size += 8; // hashval key
            size += positions.len() * std::mem::size_of::<usize>();
        }
        if let Some(ref seq) = self.raw_sequence {
            size += seq.len();
        }
        size
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SEQ: &str =
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR";

    // Exact outputs for SEQ (78 residues) at protein k=5, scaled=1: all 74 k-mers
    // hash distinctly, so both the minhash and the position map hold 74 entries.
    const SEQ_MINS: usize = 74;
    // HP-encoding of SEQ (each residue -> h/p), stored for non-protein moltypes.
    const SEQ_HP_ENCODED: &str =
        "hpphhhhpppphphhppphppphppphhhhphphhhhpphhphpppphphhpphhphphphhhphphphhpphhphpp";

    /// k=0 used to reach `add_protein`, where walking windows with
    /// `len().saturating_sub(ksize - 1)` underflows and panics. The size is now
    /// rejected at construction instead.
    #[test]
    fn test_new_rejects_zero_ksize() {
        let err = ProteinSketch::new("p", 0, 1, "protein").unwrap_err();
        assert!(
            err.to_string().contains("K-mer size must be greater than 0"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn test_new_rejects_oversized_ksize() {
        let err = ProteinSketch::new("p", 101, 1, "protein").unwrap_err();
        assert!(err.to_string().contains("K-mer size too large"), "unexpected error: {err}");
    }

    #[test]
    fn test_new_accepts_boundary_ksizes() {
        assert_eq!(ProteinSketch::new("p", 1, 1, "protein").unwrap().protein_ksize(), 1);
        assert_eq!(ProteinSketch::new("p", 100, 1, "protein").unwrap().protein_ksize(), 100);
    }

    /// k=1 is the smallest accepted size and the boundary next to the rejected 0,
    /// so it should sketch normally rather than panic.
    #[test]
    fn test_add_protein_at_ksize_one() {
        let mut s = ProteinSketch::new("p", 1, 1, "protein").unwrap();
        s.add_protein(SEQ, false).unwrap();
        // Distinct 1-mers are bounded by the alphabet, not the 78 windows. SEQ uses
        // 19 of the 20 amino acids; it contains no cysteine.
        assert_eq!(s.kmer_positions().len(), 19);
    }

    fn protein_sketch() -> ProteinSketch {
        ProteinSketch::from_protein_sequence("p1", SEQ, 5, 1, "protein20").unwrap()
    }

    #[test]
    fn test_new_empty_sketch_accessors() {
        let s = ProteinSketch::new("n", 5, 1, "hp_lehninger2").unwrap();
        assert_eq!(s.protein_ksize(), 5);
        assert_eq!(s.scaled(), 1);
        assert_eq!(s.minhash_ksize(), 15); // 5 * PROTEIN_TO_MINHASH_RATIO
                                           // `hp` normalizes to the current name for the Lehninger 2-class alphabet.
        assert_eq!(s.moltype().to_string(), "hp_lehninger2");
        assert_eq!(s.signature().minhash.mins().len(), 0);
        assert!(!s.has_efficient_data());
        assert!(s.get_efficient_data().is_none());
        assert_eq!(s.get_raw_sequence(), None);
        assert_eq!(s.get_moltype_sequence(), None);
    }

    #[test]
    fn test_from_protein_sequence_builds_minhash_positions_and_raw_seq() {
        let s = protein_sketch();
        assert_eq!(s.signature().minhash.mins().len(), SEQ_MINS);
        assert_eq!(s.kmer_positions().len(), SEQ_MINS);
        assert!(s.has_efficient_data());
        // store_sequences=true keeps the raw sequence; protein moltype stores no encoding.
        assert_eq!(s.get_raw_sequence(), Some(SEQ));
        assert_eq!(s.get_moltype_sequence(), None);
    }

    #[test]
    fn test_non_protein_moltype_stores_encoded_sequence() {
        let s = ProteinSketch::from_protein_sequence("p", SEQ, 5, 1, "hp_lehninger2").unwrap();
        assert_eq!(s.get_moltype_sequence(), Some(SEQ_HP_ENCODED));
    }

    #[test]
    fn test_remove_low_complexity_hp_kmers_toggle() {
        use crate::tests::test_fixtures::TEST_PROTEIN;

        // TEST_PROTEIN = "PLANTANDANIMALGENQMES"; the window at position 10,
        // "IMALG", is all-hydrophobic ("hhhhh") under the HP (Lehninger)
        // alphabet — the low-complexity case this targets.
        const IMALG_HASH: u64 = 8541583772724823208;

        // Default (off): low-complexity k-mers are kept, matching legacy behavior.
        let mut off = ProteinSketch::new("off", 5, 1, "hp_lehninger2").unwrap();
        off.add_protein(TEST_PROTEIN, false).unwrap();
        assert_eq!(off.kmer_positions().len(), 14);
        assert!(off.kmer_positions().contains_key(&IMALG_HASH));

        // Opted in: the all-hydrophobic "IMALG" k-mer is dropped.
        let mut on = ProteinSketch::new("on", 5, 1, "hp_lehninger2").unwrap();
        on.set_remove_low_complexity(true);
        on.add_protein(TEST_PROTEIN, false).unwrap();
        assert_eq!(on.kmer_positions().len(), 13);
        assert!(!on.kmer_positions().contains_key(&IMALG_HASH));
    }

    // Residues 26-55 of human FKBP8 (UniProt Q14318), which contains a genuine
    // 11-residue poly-glutamate (E) tract — a real low-complexity region, not
    // an invented motif.
    const FKBP8_POLY_E: &str = "VLDGVEDAEGEEEEEEEEEEEDDLSELPPL";

    #[test]
    fn test_remove_low_complexity_raw_amino_acid_homopolymer() {
        // Hash of "EEEEE" under the identity (protein) encoding.
        const POLY_E_HASH: u64 = 11331501307295692494;

        // Default (off): the raw poly-E homopolymer k-mer is kept, matching
        // legacy behavior. The 7 overlapping "EEEEE" windows (positions 10-16
        // within FKBP8_POLY_E) collapse to a single hash entry with 7 positions.
        let mut off = ProteinSketch::new("off", 5, 1, "protein20").unwrap();
        off.add_protein(FKBP8_POLY_E, false).unwrap();
        assert!(off.kmer_positions().contains_key(&POLY_E_HASH));
        assert_eq!(off.kmer_positions()[&POLY_E_HASH].len(), 7);
        assert_eq!(off.kmer_positions().len(), 20);

        // Opted in: the raw poly-E k-mer is dropped, even though this is
        // "protein" moltype (no HP encoding involved at all).
        let mut on = ProteinSketch::new("on", 5, 1, "protein20").unwrap();
        on.set_remove_low_complexity(true);
        on.add_protein(FKBP8_POLY_E, false).unwrap();
        assert!(!on.kmer_positions().contains_key(&POLY_E_HASH));
        assert_eq!(on.kmer_positions().len(), 19);
    }

    /// Custom `hp_*` alphabets take a separate encode-and-check path from the
    /// sourmash built-in `hp` (they pre-encode to uppercase H/P via their own
    /// table), so removal has to be exercised there too.
    #[test]
    fn test_remove_low_complexity_custom_hp_alphabet() {
        use crate::tests::test_fixtures::TEST_PROTEIN;

        // "IMALG" (position 10 of TEST_PROTEIN) is all-hydrophobic under the
        // Lehninger partition, which places G in the h class. It is not a raw
        // amino-acid homopolymer, so only the HP-encoded check can catch it --
        // the branch this test covers.
        let mut off = ProteinSketch::new("off", 5, 1, "hp_lehninger2").unwrap();
        off.add_protein(TEST_PROTEIN, false).unwrap();
        assert_eq!(off.kmer_positions().len(), 14);
        assert_eq!(off.low_complexity_counts(), (0, 0), "counters stay 0 when removal is off");

        let mut on = ProteinSketch::new("on", 5, 1, "hp_lehninger2").unwrap();
        on.set_remove_low_complexity(true);
        on.add_protein(TEST_PROTEIN, false).unwrap();
        assert_eq!(on.kmer_positions().len(), 13);

        // 21 residues at k=5 gives 17 windows; exactly one ("IMALG") is removed.
        assert_eq!(on.low_complexity_counts(), (17, 1));
    }

    /// A raw homopolymer and a window that only becomes a run once encoded are dropped
    /// by the same check, since a raw run always encodes to a run.
    #[test]
    fn test_remove_low_complexity_counts_raw_and_encoded_runs_together() {
        let mut on = ProteinSketch::new("on", 5, 1, "hp_lehninger2").unwrap();
        on.set_remove_low_complexity(true);
        on.add_protein(FKBP8_POLY_E, false).unwrap();

        // FKBP8_POLY_E is 30 residues -> 26 windows at k=5. Its HP (Lehninger)
        // encoding is:
        //   VLDGVEDAEGEEEEEEEEEEEDDLSELPPL
        //   hhphhpphphppppppppppppphpphhhh
        // The 11-residue E tract yields 7 fully-inside "EEEEE" windows. Two further
        // windows ("EEEED", "EEEDD") are not raw homopolymers but still encode to
        // "ppppp". 7 + 2 = 9.
        assert_eq!(on.low_complexity_counts(), (26, 9));
    }

    /// dayhoff6 is encoded by sourmash rather than by a table of ours, and used to get
    /// only the raw check, so `EEEED` and `EEEDD` (both `ccccc`: D and E share the
    /// Dayhoff acid/amide class) were kept while the same windows under an HP alphabet
    /// were dropped.
    #[test]
    fn test_remove_low_complexity_checks_encoded_runs_under_dayhoff() {
        let mut on = ProteinSketch::new("on", 5, 1, "dayhoff6").unwrap();
        on.set_remove_low_complexity(true);
        on.add_protein(FKBP8_POLY_E, false).unwrap();

        // Dayhoff encoding of FKBP8_POLY_E:
        //   VLDGVEDAEGEEEEEEEEEEEDDLSELPPL
        //   eecbeccbcbcccccccccccccebcebbe
        // The 13-residue run of c (E x 11 then DD) holds 13 - 5 + 1 = 9 windows.
        assert_eq!(on.low_complexity_counts(), (26, 9));
        assert_eq!(on.kmer_positions().len(), 17);
    }

    #[test]
    fn test_efficient_data_roundtrip() {
        let s = protein_sketch();
        let store = s.to_efficient_data(true);
        assert_eq!(store.name, "p1");
        assert_eq!(store.kmer_count(), SEQ_MINS);
        assert_eq!(store.get_raw_sequence(), Some(SEQ));
        // name(2) + 74 mins*8 + 74 abunds*8 + 74 pos entries*(8 + 1*8) + raw seq(78) = 2448.
        assert_eq!(store.estimated_size(), 2448);

        let restored =
            ProteinSketch::from_efficient_data(store, "protein20".to_string(), 5, 1).unwrap();
        assert_eq!(restored.signature().minhash.mins(), s.signature().minhash.mins());
        assert_eq!(restored.kmer_positions(), s.kmer_positions());
        assert!(restored.has_efficient_data());
    }

    #[test]
    fn test_serde_json_roundtrip() {
        let s = protein_sketch();
        let json = serde_json::to_string(&s).unwrap();
        let back: ProteinSketch = serde_json::from_str(&json).unwrap();
        assert_eq!(back.protein_ksize(), s.protein_ksize());
        assert_eq!(back.scaled(), s.scaled());
        assert_eq!(back.moltype(), s.moltype());
        assert_eq!(back.signature().minhash.mins(), s.signature().minhash.mins());
        assert_eq!(back.kmer_positions(), s.kmer_positions());
        // efficient_data is not part of the serialized form.
        assert!(!back.has_efficient_data());
    }

    #[test]
    fn test_compatibility_and_set_ops() {
        let a = protein_sketch();
        let b = protein_sketch();
        assert!(a.is_compatible(&b));
        assert!(!a.is_compatible(&ProteinSketch::new("c", 6, 1, "protein20").unwrap()));
        assert!(!a.is_compatible(&ProteinSketch::new("d", 5, 1, "hp_lehninger2").unwrap()));

        let mins = a.mins_as_set();
        assert_eq!(mins.len(), SEQ_MINS);
        // Identical sketches intersect fully.
        assert_eq!(a.intersect(&b).len(), SEQ_MINS);
    }

    #[test]
    fn test_set_and_get_efficient_data() {
        let mut s = ProteinSketch::new("n", 5, 1, "protein20").unwrap();
        assert!(!s.has_efficient_data());
        let store =
            ProteinSketchStore::new("n".into(), vec![1, 2], None, HashMap::new(), None, None);
        s.set_efficient_data(store);
        assert!(s.has_efficient_data());
        assert!(s.get_efficient_data().is_some());
    }

    #[test]
    fn test_into_signature() {
        let s = protein_sketch();
        let mins = s.signature().minhash.mins();
        assert_eq!(s.into_signature().minhash.mins(), mins);
    }

    #[test]
    fn test_store_new_getters_and_estimated_size() {
        let mut store = ProteinSketchStore::new(
            "n".into(),
            vec![1, 2, 3],
            Some(vec![1, 1, 1]),
            HashMap::new(),
            None,
            None,
        );
        assert_eq!(store.kmer_count(), 3);
        assert!(!store.has_raw_sequence_storage());
        assert!(store.get_raw_sequence().is_none());
        // set_raw_sequence is a no-op when raw storage was not allocated.
        store.set_raw_sequence("IGNORED".into());
        assert!(store.get_raw_sequence().is_none());
        store.set_encoded_sequence("ENCODED".into());
        assert_eq!(store.get_encoded_sequence(), Some("ENCODED"));
        // name "n"(1) + 3 mins*8 + 3 abunds*8, no positions/raw seq = 49.
        assert_eq!(store.estimated_size(), 49);
    }

    #[test]
    fn test_store_with_sequence_capacity_enables_raw_storage() {
        let mut store = ProteinSketchStore::with_sequence_capacity(
            "n".into(),
            vec![9],
            None,
            HashMap::new(),
            16,
        );
        assert!(store.has_raw_sequence_storage());
        store.set_raw_sequence("HELLO".into());
        assert_eq!(store.get_raw_sequence(), Some("HELLO"));

        // Zero capacity leaves raw storage unallocated.
        let store0 = ProteinSketchStore::with_sequence_capacity(
            "n".into(),
            vec![9],
            None,
            HashMap::new(),
            0,
        );
        assert!(!store0.has_raw_sequence_storage());
    }
}
