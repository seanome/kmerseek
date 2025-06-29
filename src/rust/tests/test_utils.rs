use crate::signature::ProteinSignature;

/// Utility function to print k-mer infos in a readable format for debugging
pub fn print_kmer_infos(protein_signature: &ProteinSignature) {
    println!("\n---\nmd5sum:  {}", protein_signature.signature().md5sum);
    println!("Name: {}", protein_signature.signature().name);
    println!("Len of Kmer infos: {}", protein_signature.kmer_infos().len());
    println!("Kmer Info Details:");
    println!("Hash\t\t\t{}\tOrig\tPos", protein_signature.moltype());
    println!("----------------------------------------");
    for (hash, kmer_info) in protein_signature.kmer_infos().iter() {
        for (original_kmer, positions) in &kmer_info.original_kmer_to_position {
            println!("{}\t{}\t{}\t{:?}", hash, kmer_info.encoded_kmer, original_kmer, positions);
        }
    }
    println!("----------------------------------------\n");
}
