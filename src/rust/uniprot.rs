use anyhow::{bail, Result};
use flate2::read::GzDecoder;
use quick_xml::events::Event;
use quick_xml::reader::Reader;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;

/// Represents a protein feature or domain from UniProt
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ProteinFeature {
    pub feature_type: String,     // e.g. "repeat", "region of interest", "chain"
    pub description: String,      // e.g. "HAT 1", "Interaction with CETN2"
    pub start: usize,             // 1-based position start
    pub end: usize,               // 1-based position end
    pub evidence: Option<String>, // Evidence code from UniProt
}

/// Represents a protein sequence and its features from UniProt
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct UniProtEntry {
    pub id: String,        // e.g. "SFI1_CALJA"
    pub accession: String, // e.g. "B0KWR6"
    pub sequence: String,  // The amino acid sequence
    pub features: Vec<ProteinFeature>,
}

impl Default for UniProtEntry {
    fn default() -> Self {
        Self {
            id: String::new(),
            accession: String::new(),
            sequence: String::new(),
            features: Vec::new(),
        }
    }
}

impl UniProtEntry {
    /// Parse a UniProt XML file and extract protein information
    pub fn from_xml<P: AsRef<Path>>(path: P) -> Result<Vec<UniProtEntry>> {
        let path_str = path.as_ref().to_string_lossy();
        let file = File::open(&path)?;
        let reader = BufReader::new(file);

        // If the file ends with .gz, use GzDecoder
        let reader: Box<dyn io::BufRead> = if path_str.ends_with(".gz") {
            Box::new(BufReader::new(GzDecoder::new(reader)))
        } else {
            Box::new(reader)
        };

        let mut xml_reader = Reader::from_reader(reader);
        xml_reader.config_mut().trim_text(true);

        let mut proteins = Vec::new();
        let mut buf = Vec::new();

        let mut current_protein: Option<UniProtEntry> = None;
        let mut in_sequence = false;
        let mut sequence_buf = String::new();
        let mut in_recommended_name = false;
        let mut in_full_name = false;

        loop {
            match xml_reader.read_event_into(&mut buf) {
                Ok(Event::Start(ref e)) => {
                    match e.name().as_ref() {
                        b"entry" => {
                            current_protein = Some(UniProtEntry {
                                id: String::new(),
                                accession: String::new(),
                                sequence: String::new(),
                                features: Vec::new(),
                            });
                        }
                        b"sequence" => {
                            in_sequence = true;
                        }
                        b"recommendedName" => {
                            in_recommended_name = true;
                        }
                        b"fullName" => {
                            if in_recommended_name {
                                in_full_name = true;
                            }
                        }
                        b"feature" => {
                            if let Some(ref mut protein) = current_protein {
                                let mut feature = ProteinFeature {
                                    feature_type: String::new(),
                                    description: String::new(),
                                    start: 0,
                                    end: 0,
                                    evidence: None,
                                };

                                // Parse feature attributes
                                for attr in e.attributes() {
                                    if let Ok(attr) = attr {
                                        match attr.key.as_ref() {
                                            b"type" => {
                                                feature.feature_type =
                                                    String::from_utf8_lossy(&attr.value)
                                                        .into_owned()
                                            }
                                            b"description" => {
                                                feature.description =
                                                    String::from_utf8_lossy(&attr.value)
                                                        .into_owned()
                                            }
                                            _ => {}
                                        }
                                    }
                                }

                                protein.features.push(feature);
                            }
                        }
                        b"accession" => {
                            if let Some(ref mut protein) = current_protein {
                                // Parse accession
                                if let Ok(Event::Text(e)) = xml_reader.read_event_into(&mut buf) {
                                    protein.accession = e.unescape()?.into_owned();
                                }
                            }
                        }
                        _ => {}
                    }
                }
                Ok(Event::Text(e)) => {
                    if in_sequence {
                        sequence_buf.push_str(&e.unescape()?.into_owned());
                    } else if in_full_name {
                        if let Some(ref mut protein) = current_protein {
                            protein.id = e.unescape()?.into_owned();
                        }
                    }
                }
                Ok(Event::End(ref e)) => match e.name().as_ref() {
                    b"entry" => {
                        if let Some(mut protein) = current_protein.take() {
                            protein.sequence = sequence_buf.clone();
                            proteins.push(protein);
                            sequence_buf.clear();
                        }
                    }
                    b"sequence" => {
                        in_sequence = false;
                    }
                    b"recommendedName" => {
                        in_recommended_name = false;
                    }
                    b"fullName" => {
                        in_full_name = false;
                    }
                    _ => {}
                },
                Ok(Event::Eof) => break,
                Err(e) => bail!("Error parsing XML: {}", e),
                _ => {}
            }
        }

        Ok(proteins)
    }
}
