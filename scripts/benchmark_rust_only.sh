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
        # Test different alphabets
        for alphabet in protein20 hp_lehninger2 dayhoff6; do
            echo "  k=$ksize, alphabet=$alphabet:"
            
            # Time the execution
            START_TIME=$(date +%s.%N)
            
            kmerseek index \
                --input "$test_file" \
                --output "$TEMP_DIR/output_${ksize}_${alphabet}.db" \
                --ksize "$ksize" \
                --alphabet "$alphabet"
            
            END_TIME=$(date +%s.%N)
            EXECUTION_TIME=$(echo "$END_TIME - $START_TIME" | bc)
            
            # Get output size
            OUTPUT_SIZE=$(du -h "$TEMP_DIR/output_${ksize}_${alphabet}.db" | cut -f1)
            
            echo "    Time: ${EXECUTION_TIME}s, Size: ${OUTPUT_SIZE}"
        done
    done
    echo
done

echo "=== Memory Usage Test ==="
echo

# Test memory usage for the larger file
LARGE_FILE="tests/testdata/fasta/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz"

echo "Memory usage for large file (k=10, alphabet=hp_lehninger2):"
echo "Rust CLI:"
/usr/bin/time -l kmerseek index \
    --input "$LARGE_FILE" \
    --output "$TEMP_DIR/memory_test.db" \
    --ksize 10 \
    --alphabet hp_lehninger2 2>&1 | grep "maximum resident set size"

echo
echo "=== Benchmark Complete ===" 