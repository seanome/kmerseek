use crate::sketch::ProteinSketch;

/// Utility function to print k-mer positions for debugging
pub fn print_kmer_positions(protein_signature: &ProteinSketch) {
    println!("\n---\nmd5sum:  {}", protein_signature.signature().md5sum);
    println!("Name: {}", protein_signature.signature().name);
    println!("Num distinct hashes: {}", protein_signature.kmer_positions().len());
    println!("Hash\t\t\t\tPositions");
    println!("----------------------------------------");
    let mut entries: Vec<_> = protein_signature.kmer_positions().iter().collect();
    entries.sort_by_key(|(h, _)| *h);
    for (hash, positions) in entries {
        println!("{}\t{:?}", hash, positions);
    }
    println!("----------------------------------------\n");
}
