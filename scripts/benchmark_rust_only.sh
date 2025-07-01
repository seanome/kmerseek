#!/bin/bash

# Rust CLI Performance Benchmark
# Usage: ./scripts/benchmark_rust_only.sh

set -e

echo "=== Rust CLI Performance Benchmark ==="
echo

# Test files
TEST_FILES=(
    "tests/testdata/fasta/ced9.fasta:small"
    "tests/testdata/fasta/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz:medium"
)

# Create temporary directory
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

for test_file_info in "${TEST_FILES[@]}"; do
    IFS=':' read -r test_file size_label <<< "$test_file_info"
    
    echo "Testing $size_label file: $(basename "$test_file")"
    echo "----------------------------------------"
    
    # Test different k-mer sizes
    for ksize in 5 10 15; do
        # Test different encodings
        for encoding in protein hp dayhoff; do
            echo "  k=$ksize, encoding=$encoding:"
            
            # Time the execution
            START_TIME=$(date +%s.%N)
            
            kmerseek-rust index \
                --input "$test_file" \
                --output "$TEMP_DIR/output_${ksize}_${encoding}.db" \
                --ksize "$ksize" \
                --encoding "$encoding"
            
            END_TIME=$(date +%s.%N)
            EXECUTION_TIME=$(echo "$END_TIME - $START_TIME" | bc)
            
            # Get output size
            OUTPUT_SIZE=$(du -h "$TEMP_DIR/output_${ksize}_${encoding}.db" | cut -f1)
            
            echo "    Time: ${EXECUTION_TIME}s, Size: ${OUTPUT_SIZE}"
        done
    done
    echo
done

echo "=== Memory Usage Test ==="
echo

# Test memory usage for the larger file
LARGE_FILE="tests/testdata/fasta/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz"

echo "Memory usage for large file (k=10, encoding=hp):"
echo "Rust CLI:"
/usr/bin/time -l kmerseek-rust index \
    --input "$LARGE_FILE" \
    --output "$TEMP_DIR/memory_test.db" \
    --ksize 10 \
    --encoding hp 2>&1 | grep "maximum resident set size"

echo
echo "=== Benchmark Complete ===" 