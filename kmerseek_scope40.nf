#!/usr/bin/env nextflow

/*
 * Nextflow pipeline to run kmerseek on SCOPe40 database with various encodings and k-sizes
 *
 * Usage:
 *   nextflow run kmerseek_scope40.nf
 */

// FASTA files to process
// fasta_gd: Standard amino acid sequences
// fasta_tea: Same sequences with "The Embedded Alphabet" encoding
params.fasta_gd = "${System.getProperty('user.home')}/data/scope/astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.fa"
params.fasta_tea = "${System.getProperty('user.home')}/data/scope/astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.tea.fa"
params.outdir = "${System.getProperty('user.home')}/data/scope/results-2025-12-25-average_kmer_rarity"
params.kmerseek = "${System.getProperty('user.home')}/code/kmerseek/target/release/kmerseek-rust"

// Define k-size ranges for each encoding
// hp: 15-30
// protein: 5-10
// dayhoff: 10-15 (reduced range to avoid memory issues with smaller k-sizes)

process buildRelease {
    output:
    path 'kmerseek-rust'

    script:
    def kmerseek_dir = "${System.getProperty('user.home')}/code/kmerseek"
    """
    WORK_DIR=\$PWD
    cd ${kmerseek_dir}
    cargo build --release
    cp ${kmerseek_dir}/target/release/kmerseek-rust \$WORK_DIR/kmerseek-rust
    """
}

process indexDatabase {
    tag "${encoding}_k${ksize}_scaled${scaled}_${fasta.name.replaceAll('\\.fa$', '').replaceAll('\\.fasta$', '')}"
    publishDir params.outdir, mode: 'copy', pattern: '*.rocksdb/**'
    publishDir params.outdir, mode: 'copy', pattern: '*.index.log'

    input:
    path kmerseek
    tuple path(fasta), val(encoding), val(ksize), val(scaled)

    output:
    tuple val(encoding), val(ksize), val(scaled), path("${fasta}.${encoding}.k${ksize}.scaled${scaled}.kmerseek.rocksdb")
    path "${fasta}.${encoding}.k${ksize}.scaled${scaled}.kmerseek.index.log"

    script:
    def log_file = "${fasta}.${encoding}.k${ksize}.scaled${scaled}.kmerseek.index.log"
    """
    echo "=== Indexing: ${encoding} k=${ksize} scaled=${scaled} ===" | tee ${log_file}
    echo "Start time: \$(date '+%Y-%m-%d %H:%M:%S')" | tee -a ${log_file}
    echo "" | tee -a ${log_file}

    /usr/bin/time -l ${kmerseek} index \\
        --encoding ${encoding} \\
        --ksize ${ksize} \\
        --scaled ${scaled} \\
        --input ${fasta} \\
        2>&1 | tee -a ${log_file}

    echo "" | tee -a ${log_file}
    echo "End time: \$(date '+%Y-%m-%d %H:%M:%S')" | tee -a ${log_file}
    """
}

process searchAllVsAll {
    tag "${encoding}_k${ksize}_scaled${scaled}_${index.name.replaceAll('\\.rocksdb$', '')}"
    publishDir params.outdir, mode: 'copy', pattern: '*.csv'
    publishDir params.outdir, mode: 'copy', pattern: '*.search.log'

    input:
    path kmerseek
    tuple val(encoding), val(ksize), val(scaled), path(index)

    output:
    tuple val(encoding), val(ksize), val(scaled), path("${index.name.replaceAll('\\.rocksdb$', '')}.results.csv")
    path "${index.name.replaceAll('\\.rocksdb$', '')}.search.log"

    script:
    def basename = index.name.replaceAll("\\.rocksdb\$", '')
    def output_csv = "${basename}.results.csv"
    def log_file = "${basename}.search.log"
    """
    echo "=== Searching: ${encoding} k=${ksize} ===" | tee ${log_file}
    echo "Start time: \$(date '+%Y-%m-%d %H:%M:%S')" | tee -a ${log_file}
    echo "" | tee -a ${log_file}

    /usr/bin/time -l ${kmerseek} search \\
        --encoding ${encoding} \\
        --ksize ${ksize} \\
        --query-is-index \\
        --query ${index} \\
        --target ${index} \\
        > ${output_csv} 2>> ${log_file}

    echo "" | tee -a ${log_file}
    echo "End time: \$(date '+%Y-%m-%d %H:%M:%S')" | tee -a ${log_file}
    echo "Results: \$(wc -l < ${output_csv}) rows" | tee -a ${log_file}
    """
}

workflow {
    // Build kmerseek in release mode
    kmerseek_bin = buildRelease()

    // Create k-size channels
    // HP k=15-45 (re-running with rebuilt kmerseek to get new columns)
    hp_ksizes = Channel.of(15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45)
    protein_ksizes = Channel.of(5, 6, 7, 8, 9, 10)
    dayhoff_ksizes = Channel.of(10, 11, 12, 13, 14, 15)

    // Create FASTA file channels
    fasta_gd_ch = Channel.fromPath(params.fasta_gd)
    fasta_tea_ch = Channel.fromPath(params.fasta_tea)

    // Create encoding-ksize-scaled combinations for standard sequences (scaled=1 for all)
    hp_params = hp_ksizes.map { k -> tuple('hp', k, 1) }
    protein_params_gd = protein_ksizes.map { k -> tuple('protein', k, 1) }
    dayhoff_params = dayhoff_ksizes.map { k -> tuple('dayhoff', k, 1) }
    all_params_gd = hp_params
    //.mix(protein_params_gd, dayhoff_params)

    // Create encoding-ksize-scaled combinations for TEA-encoded sequences (protein k=20 with scaled=10)
    // TEA = "The Embedded Alphabet" - same sequences, different character encoding
    // Memory usage with scaled=1: k=8: 107 GB, k=10: 92 GB, k=15: 68 GB, k=20: 64 GB, k=30: 39 GB (all crashed)
    // Memory usage with scaled=10: k=10: 94 GB (crashed - somehow MORE than scaled=1!)
    // Trying k=20 with scaled=10 (k=20 scaled=1 was 64GB, so scaled=10 should be ~6GB)
    protein_ksizes_tea = Channel.of(20)
    protein_params_tea = protein_ksizes_tea.map { k -> tuple('protein', k, 10) }

    // Combine FASTA files with their respective parameters
    fasta_params_gd = fasta_gd_ch.combine(all_params_gd)
    fasta_params_tea = fasta_tea_ch.combine(protein_params_tea)

    // Running HP k=15-45 on standard sequences (re-running with rebuilt kmerseek for new columns)
    all_fasta_params = fasta_params_gd  // TEA later: fasta_params_gd.mix(fasta_params_tea)

    // Index the database for each encoding/ksize combination
    indexed = indexDatabase(kmerseek_bin, all_fasta_params)

    // Extract just the index tuple (encoding, ksize, scaled, path) from the output
    // indexDatabase outputs both the index tuple and the log file
    index_only = indexed[0]

    // Run all-vs-all search on each index
    // searchAllVsAll outputs both the csv tuple and the log file
    search_outputs = searchAllVsAll(kmerseek_bin, index_only)

    // Extract just the CSV results tuple
    results = search_outputs[0]

    // Print summary
    results.subscribe { encoding, ksize, scaled, csv ->
        println("Completed: ${encoding} k=${ksize} scaled=${scaled} -> ${csv}")
    }
}
