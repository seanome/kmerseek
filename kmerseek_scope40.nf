#!/usr/bin/env nextflow

/*
 * Nextflow pipeline to run kmerseek on SCOPe40 database with various encodings and k-sizes
 *
 * Usage:
 *   nextflow run kmerseek_scope40.nf
 */

params.fasta = "${System.getProperty('user.home')}/data/scope/astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.fa"
params.outdir = "${System.getProperty('user.home')}/data/scope/results-2025-12-25-average_kmer_rarity"
params.kmerseek = "${System.getProperty('user.home')}/code/kmerseek/target/release/kmerseek-rust"

// Define k-size ranges for each encoding
// hp: 16-20
// protein: 5-10
// dayhoff: 10-15 (reduced range to avoid memory issues with smaller k-sizes)
hp_ksizes = Channel.of(15, 16, 17, 18, 19, 20)
protein_ksizes = Channel.of(5, 6, 7, 8, 9, 10)
dayhoff_ksizes = Channel.of(10, 11, 12, 13, 14, 15)

// Create encoding-ksize combinations
hp_params = hp_ksizes.map { k -> tuple('hp', k) }
protein_params = protein_ksizes.map { k -> tuple('protein', k) }
dayhoff_params = dayhoff_ksizes.map { k -> tuple('dayhoff', k) }

// Combine all parameters
all_params = hp_params.mix(protein_params, dayhoff_params)

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
    tag "${encoding}_k${ksize}"
    publishDir params.outdir, mode: 'copy', pattern: '*.rocksdb/**'
    publishDir params.outdir, mode: 'copy', pattern: '*.index.log'

    input:
    path kmerseek
    tuple path(fasta), val(encoding), val(ksize)

    output:
    tuple val(encoding), val(ksize), path("${fasta}.${encoding}.k${ksize}.scaled1.kmerseek.rocksdb")
    path "${fasta}.${encoding}.k${ksize}.scaled1.kmerseek.index.log"

    script:
    def log_file = "${fasta}.${encoding}.k${ksize}.scaled1.kmerseek.index.log"
    """
    echo "=== Indexing: ${encoding} k=${ksize} ===" | tee ${log_file}
    echo "Start time: \$(date '+%Y-%m-%d %H:%M:%S')" | tee -a ${log_file}
    echo "" | tee -a ${log_file}

    /usr/bin/time -l ${kmerseek} index \\
        --encoding ${encoding} \\
        --ksize ${ksize} \\
        --input ${fasta} \\
        2>&1 | tee -a ${log_file}

    echo "" | tee -a ${log_file}
    echo "End time: \$(date '+%Y-%m-%d %H:%M:%S')" | tee -a ${log_file}
    """
}

process searchAllVsAll {
    tag "${encoding}_k${ksize}"
    publishDir params.outdir, mode: 'copy', pattern: '*.csv'
    publishDir params.outdir, mode: 'copy', pattern: '*.search.log'

    input:
    path kmerseek
    tuple val(encoding), val(ksize), path(index)

    output:
    tuple val(encoding), val(ksize), path("${index.name.replaceAll('\\.rocksdb$', '')}.results.csv")
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

    // Create input channel with fasta file
    fasta_ch = Channel.fromPath(params.fasta)

    // Combine fasta with all parameters to create all combinations
    fasta_params = fasta_ch.combine(all_params)

    // Index the database for each encoding/ksize combination
    indexed = indexDatabase(kmerseek_bin, fasta_params)

    // Extract just the index tuple (encoding, ksize, path) from the output
    // indexDatabase outputs both the index tuple and the log file
    index_only = indexed[0]

    // Run all-vs-all search on each index
    // searchAllVsAll outputs both the csv tuple and the log file
    search_outputs = searchAllVsAll(kmerseek_bin, index_only)

    // Extract just the CSV results tuple
    results = search_outputs[0]

    // Print summary
    results.subscribe { encoding, ksize, csv ->
        println "Completed: ${encoding} k=${ksize} -> ${csv}"
    }
}
