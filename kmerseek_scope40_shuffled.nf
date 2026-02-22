#!/usr/bin/env nextflow

/*
 * Nextflow pipeline to search shuffled SCOPe40 sequences against original SCOPe40 database
 * This tests the false positive rate of kmerseek by searching 2-mer-preserving shuffled sequences
 *
 * Usage:
 *   nextflow run kmerseek_scope40_shuffled.nf
 */

// Import processes from the main pipeline
include { buildRelease } from './kmerseek_scope40.nf'
include { indexDatabase } from './kmerseek_scope40.nf'

// FASTA files
params.fasta_shuffled = "${System.getProperty('user.home')}/data/scope/astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.shuffled_2mer_x3.fa"
// Pre-existing indices for original database (from previous pipeline run)
params.original_indices = "${System.getProperty('user.home')}/code/kmerseek/work/**/astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.fa.hp.k*.scaled1.kmerseek.rocksdb"
params.outdir = "${System.getProperty('user.home')}/data/scope/results-shuffled-2025-02-07"
params.kmerseek = "${System.getProperty('user.home')}/code/kmerseek/target/release/kmerseek-rust"

process searchShuffledVsOriginal {
    tag "${encoding}_k${ksize}_scaled${scaled}"
    publishDir params.outdir, mode: 'copy', pattern: '*.csv'
    publishDir params.outdir, mode: 'copy', pattern: '*.search.log'

    input:
    path kmerseek
    tuple val(encoding), val(ksize), val(scaled), path(query_index), path(target_index)

    output:
    tuple val(encoding), val(ksize), val(scaled), path("shuffled_vs_original.${encoding}.k${ksize}.scaled${scaled}.results.csv")
    path "shuffled_vs_original.${encoding}.k${ksize}.scaled${scaled}.search.log"

    script:
    def output_csv = "shuffled_vs_original.${encoding}.k${ksize}.scaled${scaled}.results.csv"
    def log_file = "shuffled_vs_original.${encoding}.k${ksize}.scaled${scaled}.search.log"
    """
    echo "=== Searching shuffled vs original: ${encoding} k=${ksize} scaled=${scaled} ===" | tee ${log_file}
    echo "Query: ${query_index}" | tee -a ${log_file}
    echo "Target: ${target_index}" | tee -a ${log_file}
    echo "Start time: \$(date '+%Y-%m-%d %H:%M:%S')" | tee -a ${log_file}
    echo "" | tee -a ${log_file}

    /usr/bin/time -l ${kmerseek} search \\
        --encoding ${encoding} \\
        --ksize ${ksize} \\
        --query-is-index \\
        --query ${query_index} \\
        --target ${target_index} \\
        > ${output_csv} 2>> ${log_file}

    echo "" | tee -a ${log_file}
    echo "End time: \$(date '+%Y-%m-%d %H:%M:%S')" | tee -a ${log_file}
    echo "Results: \$(wc -l < ${output_csv}) rows" | tee -a ${log_file}
    """
}

workflow {
    // Build kmerseek in release mode
    kmerseek_bin = buildRelease()

    // Load pre-existing original indices and extract encoding, ksize, scaled from filename
    // Filename pattern: astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.fa.hp.k15.scaled1.kmerseek.rocksdb
    // Use .unique() to deduplicate indices from multiple Nextflow runs
    original_indices = channel
        .fromPath(params.original_indices, type: 'dir')
        .map { path ->
            def name = path.name
            def matcher = name =~ /\.(\w+)\.k(\d+)\.scaled(\d+)\.kmerseek\.rocksdb$/
            if (matcher) {
                def encoding = matcher[0][1]
                def ksize = matcher[0][2] as Integer
                def scaled = matcher[0][3] as Integer
                return tuple(encoding, ksize, scaled, path)
            }
            return null
        }
        .filter { item -> item != null }
        .unique { encoding, ksize, scaled, _path -> tuple(encoding, ksize, scaled) }
        .filter { _encoding, ksize, _scaled, _path -> ksize >= 20 }  // Skip k<20 due to memory limits

    // Get the ksizes that have pre-existing indices
    // Create matching parameters for shuffled indexing
    shuffled_params = original_indices
        .map { encoding, ksize, scaled, _path -> tuple(encoding, ksize, scaled) }

    // Create FASTA file channel for shuffled
    fasta_shuffled_ch = channel.fromPath(params.fasta_shuffled)

    // Combine shuffled FASTA with parameters from existing indices
    fasta_params_shuffled = fasta_shuffled_ch
        .combine(shuffled_params)
        .map { fasta, encoding, ksize, scaled -> tuple(fasta, encoding, ksize, scaled) }

    // Index shuffled database
    indexed_shuffled = indexDatabase(kmerseek_bin, fasta_params_shuffled)

    // Extract just the index tuples (encoding, ksize, scaled, path)
    shuffled_indices = indexed_shuffled[0]

    // Join shuffled (query) with original (target) by encoding, ksize, scaled
    // Result: tuple(encoding, ksize, scaled, shuffled_path, original_path)
    search_pairs = shuffled_indices
        .join(original_indices, by: [0, 1, 2])
        .map { encoding, ksize, scaled, shuffled_path, original_path ->
            tuple(encoding, ksize, scaled, shuffled_path, original_path)
        }

    // Search shuffled vs original
    search_outputs = searchShuffledVsOriginal(kmerseek_bin, search_pairs)

    // Extract just the CSV results tuple
    results = search_outputs[0]

    // Print summary
    results.subscribe { encoding, ksize, scaled, csv ->
        println("Completed: ${encoding} k=${ksize} scaled=${scaled} -> ${csv}")
    }
}
