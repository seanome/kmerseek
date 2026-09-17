//! One query sequence against one target sequence: every k-mer the two share, where it
//! sits in each, and the matched regions those k-mers chain into.
//!
//! This is the data behind `kmerseek pair` and `scripts/visualize_pair.py`. Search reports
//! only the matched regions; this keeps the individual shared k-mers too, so a plot can show
//! the ones that fall on a region's diagonal next to the ones scattered elsewhere.

use std::collections::HashSet;
use std::path::Path;

use needletail::parse_fastx_file;
use serde::Serialize;

use crate::aminoacid::AminoAcidAmbiguity;
use crate::errors::{IndexError, IndexResult};
use crate::search::{find_matched_regions, MatchedRegion};
use crate::sketch::ProteinSketch;

/// One sequence of the pair, with its reduced-alphabet encoding.
#[derive(Debug, Clone, Serialize)]
pub struct PairSequence {
    /// The FASTA header after `>`, the same name search reports.
    pub name: String,
    pub sequence: String,
    /// The sequence in the reduced alphabet. Equal to `sequence` for `protein20`.
    pub encoded: String,
}

/// One k-mer both sequences contain. Positions are 0-based k-mer starts, as in the search
/// CSV. A k-mer present more than once in either sequence appears once per position pair.
#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct SharedKmer {
    pub query_pos: usize,
    pub target_pos: usize,
    /// The k-mer in the reduced alphabet, which is what the two sequences agree on.
    pub kmer: String,
    /// The residues under the k-mer in each sequence, which may differ.
    pub query_kmer: String,
    pub target_kmer: String,
}

/// A run of shared k-mers consecutive in both sequences. Half-open, 0-based, in residues.
#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct PairRegion {
    pub query_start: u32,
    pub query_end: u32,
    pub target_start: u32,
    pub target_end: u32,
    pub length: u32,
}

impl From<&MatchedRegion> for PairRegion {
    fn from(region: &MatchedRegion) -> Self {
        Self {
            query_start: region.start,
            query_end: region.end,
            target_start: region.target_start,
            target_end: region.target_end,
            length: region.length,
        }
    }
}

/// Everything `kmerseek pair` writes.
#[derive(Debug, Clone, Serialize)]
pub struct PairReport {
    pub ksize: u32,
    pub moltype: String,
    pub query: PairSequence,
    pub target: PairSequence,
    /// Sorted by query position, then target position.
    pub shared_kmers: Vec<SharedKmer>,
    /// Longest first.
    pub regions: Vec<PairRegion>,
}

impl PairReport {
    /// Pretty-printed JSON, the format `scripts/visualize_pair.py` reads.
    pub fn to_json(&self) -> IndexResult<String> {
        serde_json::to_string_pretty(self).map_err(|e| IndexError::ParseError(e.to_string()))
    }
}

/// A named sequence read from a FASTA file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    pub name: String,
    pub sequence: String,
}

/// Read one record from a FASTA file: the first, or the one whose header is `name` or
/// starts with `name` followed by a space (so `sp|P10415|BCL2_HUMAN` finds the UniProt
/// header without the description).
pub fn read_record(path: &Path, name: Option<&str>) -> IndexResult<FastaRecord> {
    let mut reader = parse_fastx_file(path).map_err(|e| IndexError::ParseError(e.to_string()))?;
    while let Some(record) = reader.next() {
        let record = record.map_err(|e| IndexError::ParseError(e.to_string()))?;
        let id = std::str::from_utf8(record.id())?.to_string();
        if name.is_none_or(|wanted| header_matches(&id, wanted)) {
            let sequence = std::str::from_utf8(&record.seq())?.to_string();
            return Ok(FastaRecord { name: id, sequence });
        }
    }
    Err(IndexError::FastaParsing(match name {
        Some(wanted) => format!("no record named '{}' in {}", wanted, path.display()),
        None => format!("no records in {}", path.display()),
    }))
}

fn header_matches(header: &str, wanted: &str) -> bool {
    header == wanted || header.split_whitespace().next() == Some(wanted)
}

/// Compare two sequences at one k-mer size and alphabet. Every k-mer is kept (scaled 1).
pub fn compare_pair(
    query: &FastaRecord,
    target: &FastaRecord,
    ksize: u32,
    moltype: &str,
) -> IndexResult<PairReport> {
    let query_sketch = sketch(query, ksize, moltype)?;
    let target_sketch = sketch(target, ksize, moltype)?;
    let intersection = query_sketch.intersect(&target_sketch);
    let mut regions: Vec<PairRegion> =
        find_matched_regions(&query_sketch, &target_sketch, &intersection)
            .iter()
            .map(PairRegion::from)
            .collect();
    regions.sort_by_key(|r| (std::cmp::Reverse(r.length), r.query_start, r.target_start));
    let query = pair_sequence(&query_sketch);
    let target = pair_sequence(&target_sketch);
    let shared_kmers = shared_kmers(&query_sketch, &target_sketch, &query, &target, &intersection);
    Ok(PairReport {
        ksize,
        moltype: query_sketch.moltype().to_string(),
        query,
        target,
        shared_kmers,
        regions,
    })
}

fn sketch(record: &FastaRecord, ksize: u32, moltype: &str) -> IndexResult<ProteinSketch> {
    let sequence = AminoAcidAmbiguity::new().validate_and_resolve(&record.sequence, moltype)?;
    let mut sketch = ProteinSketch::new(&record.name, ksize, 1, moltype)?;
    sketch.add_protein(&sequence, true)?;
    Ok(sketch)
}

fn pair_sequence(sketch: &ProteinSketch) -> PairSequence {
    let sequence = sketch.get_raw_sequence().unwrap_or_default().to_string();
    let encoded = sketch.get_moltype_sequence().unwrap_or(&sequence).to_string();
    PairSequence { name: sketch.signature().name.clone(), sequence, encoded }
}

/// Every (query position, target position) pair for every shared k-mer hash.
fn shared_kmers(
    query_sketch: &ProteinSketch,
    target_sketch: &ProteinSketch,
    query: &PairSequence,
    target: &PairSequence,
    intersection: &HashSet<u64>,
) -> Vec<SharedKmer> {
    let ksize = query_sketch.protein_ksize() as usize;
    let mut shared = Vec::new();
    for hashval in intersection {
        let (Some(query_positions), Some(target_positions)) = (
            query_sketch.kmer_positions().get(hashval),
            target_sketch.kmer_positions().get(hashval),
        ) else {
            continue;
        };
        for &query_pos in query_positions {
            for &target_pos in target_positions {
                shared.push(SharedKmer {
                    query_pos,
                    target_pos,
                    kmer: query.encoded[query_pos..query_pos + ksize].to_string(),
                    query_kmer: query.sequence[query_pos..query_pos + ksize].to_string(),
                    target_kmer: target.sequence[target_pos..target_pos + ksize].to_string(),
                });
            }
        }
    }
    shared.sort_by_key(|s| (s.query_pos, s.target_pos));
    shared.dedup();
    shared
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::test_fixtures::{TEST_BLC2_FASTA, TEST_CED9_FASTA};

    fn bcl2_vs_ced9(ksize: u32, moltype: &str) -> PairReport {
        let query = read_record(Path::new(TEST_BLC2_FASTA), None).unwrap();
        let target = read_record(Path::new(TEST_CED9_FASTA), None).unwrap();
        compare_pair(&query, &target, ksize, moltype).unwrap()
    }

    #[test]
    fn read_record_takes_the_first_record_by_default() {
        let record = read_record(Path::new(TEST_BLC2_FASTA), None).unwrap();
        assert_eq!(
            record.name,
            "sp|P10415|BCL2_HUMAN Apoptosis regulator Bcl-2 OS=Homo sapiens OX=9606 GN=BCL2 PE=1 SV=2"
        );
        assert_eq!(record.sequence.len(), 239);
        assert!(record.sequence.starts_with("MAHAGRTGYDNREIVMKYIHYKLSQRGYEWDAGD"));
    }

    #[test]
    fn read_record_matches_the_accession_token() {
        let record = read_record(Path::new(TEST_CED9_FASTA), Some("sp|P41958|CED9_CAEEL")).unwrap();
        assert_eq!(record.sequence.len(), 280);
    }

    #[test]
    fn read_record_reports_a_missing_name() {
        let err = read_record(Path::new(TEST_CED9_FASTA), Some("not_here")).unwrap_err();
        assert!(err.to_string().contains("no record named 'not_here'"), "{err}");
    }

    #[test]
    fn bcl2_ced9_bh1_region_at_hp_k12() {
        let report = bcl2_vs_ced9(12, "hp");
        assert_eq!(report.moltype, "hp_lehninger2");
        assert_eq!(report.query.encoded.len(), 239);
        // The BH1 motif: BCL-2 NWGR at 143 (1-based) lines up with CED-9 SYGR at 167.
        let bh1 = &report.regions[0];
        assert_eq!(
            *bh1,
            PairRegion {
                query_start: 138,
                query_end: 157,
                target_start: 162,
                target_end: 181,
                length: 19,
            }
        );
        assert_eq!(&report.query.sequence[138..157], "RDGVNWGRIVAFFEFGGVM");
        assert_eq!(&report.target.sequence[162..181], "QCPMSYGRLIGLISFGGFV");
        assert_eq!(&report.query.encoded[138..157], "pphhphhphhhhhphhhhh");
        assert_eq!(&report.target.encoded[162..181], "pphhphhphhhhhphhhhh");
    }

    #[test]
    fn bcl2_ced9_shared_kmers_are_sorted_and_cover_the_bh1_run() {
        let report = bcl2_vs_ced9(12, "hp");
        let on_bh1: Vec<&SharedKmer> = report
            .shared_kmers
            .iter()
            .filter(|s| (138..146).contains(&s.query_pos) && s.target_pos == s.query_pos + 24)
            .collect();
        // A 19-residue run holds 19 - 12 + 1 = 8 consecutive 12-mers.
        assert_eq!(on_bh1.len(), 8);
        assert_eq!(
            on_bh1[0],
            &SharedKmer {
                query_pos: 138,
                target_pos: 162,
                kmer: "pphhphhphhhh".to_string(),
                query_kmer: "RDGVNWGRIVAF".to_string(),
                target_kmer: "QCPMSYGRLIGL".to_string(),
            }
        );
        let positions: Vec<(usize, usize)> =
            report.shared_kmers.iter().map(|s| (s.query_pos, s.target_pos)).collect();
        let mut sorted = positions.clone();
        sorted.sort();
        sorted.dedup();
        assert_eq!(positions, sorted);
    }

    #[test]
    fn bcl2_ced9_shares_27_hp_12mers_and_none_at_k4_protein20() {
        assert_eq!(bcl2_vs_ced9(12, "hp").shared_kmers.len(), 27);
        assert_eq!(bcl2_vs_ced9(4, "protein20").shared_kmers.len(), 0);
    }

    #[test]
    fn protein20_encodes_to_itself_and_repeats_give_one_row_per_position_pair() {
        let report = bcl2_vs_ced9(3, "protein20");
        assert_eq!(report.moltype, "protein20");
        assert_eq!(report.query.encoded, report.query.sequence);
        // APG occurs twice in BCL-2 (positions 44 and 76) and once in CED-9 (101).
        let apg: Vec<(usize, usize)> = report
            .shared_kmers
            .iter()
            .filter(|s| s.kmer == "APG")
            .map(|s| (s.query_pos, s.target_pos))
            .collect();
        assert_eq!(apg, vec![(44, 101), (76, 101)]);
        assert_eq!(report.shared_kmers.len(), 10);
        assert_eq!(report.shared_kmers[0].query_kmer, "APG");
        assert_eq!(report.shared_kmers[0].target_kmer, "APG");
    }
}
