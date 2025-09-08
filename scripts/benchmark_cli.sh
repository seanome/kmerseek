#!/bin/bash

# Quick benchmark script to compare Rust vs Python CLI performance
# Usage: ./scripts/benchmark_cli.sh

set -e

echo "=== CLI Performance Benchmark ==="
echo "Comparing kmerseek-rust vs kmerseek (Python)"
echo

# Test file
TEST_FILE="tests/testdata/fasta/ced9.fasta"
TEST_FILE_GZ="tests/testdata/fasta/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz"
TEST_FILE_LARGE="tests/testdata/fasta/uniprotkb_protein_name_Uncharacterized_2025_04_15.fasta.gz"

# Create temporary directory
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

echo "Test 1: Small file (ced9.fasta)"
echo "--------------------------------"

# Rust version
echo "Rust CLI:"
RUST_START=$(date +%s.%N)
kmerseek-rust index --input "$TEST_FILE" --output "$TEMP_DIR/rust_output.db" --ksize 10 --encoding protein
RUST_END=$(date +%s.%N)
RUST_TIME=$(echo "$RUST_END - $RUST_START" | bc)
RUST_SIZE=$(du -h "$TEMP_DIR/rust_output.db" | cut -f1)

# Python version
echo "Python CLI:"
PYTHON_START=$(date +%s.%N)
kmerseek index "$TEST_FILE" --extract-kmers --moltype protein --ksize 10 --scaled 1 --force
PYTHON_END=$(date +%s.%N)
PYTHON_TIME=$(echo "$PYTHON_END - $PYTHON_START" | bc)
PYTHON_SIZE=$(du -h "$TEST_FILE.protein.k10.scaled1.sig.zip" | cut -f1)

echo
echo "Results for small file:"
echo "  Rust:   ${RUST_TIME}s, ${RUST_SIZE}"
echo "  Python: ${PYTHON_TIME}s, ${PYTHON_SIZE}"
echo "  Speedup: $(echo "scale=2; $PYTHON_TIME / $RUST_TIME" | bc)x"
echo

echo "Test 2: Larger gzipped file"
echo "----------------------------"

# Clean up previous outputs
rm -f "$TEST_FILE_GZ.protein.k10.scaled1.sig.zip"
rm -rf "$TEMP_DIR/rust_output2.db"

# Rust version
echo "Rust CLI:"
RUST_START=$(date +%s.%N)
kmerseek-rust index --input "$TEST_FILE_GZ" --output "$TEMP_DIR/rust_output2.db" --ksize 10 --encoding protein
RUST_END=$(date +%s.%N)
RUST_TIME=$(echo "$RUST_END - $RUST_START" | bc)
RUST_SIZE=$(du -h "$TEMP_DIR/rust_output2.db" | cut -f1)

# Python version
echo "Python CLI:"
PYTHON_START=$(date +%s.%N)
kmerseek index "$TEST_FILE_GZ" --extract-kmers --moltype protein --ksize 10 --scaled 1 --force
PYTHON_END=$(date +%s.%N)
PYTHON_TIME=$(echo "$PYTHON_END - $PYTHON_START" | bc)
PYTHON_SIZE=$(du -h "$TEST_FILE_GZ.protein.k10.scaled1.sig.zip" | cut -f1)

echo
echo "Results for larger file:"
echo "  Rust:   ${RUST_TIME}s, ${RUST_SIZE}"
echo "  Python: ${PYTHON_TIME}s, ${PYTHON_SIZE}"
echo "  Speedup: $(echo "scale=2; $PYTHON_TIME / $RUST_TIME" | bc)x"
echo

echo "Test 3: Large uncharacterized protein dataset (~2800 sequences)"
echo "---------------------------------------------------------------"

# Clean up previous outputs
rm -f "$TEST_FILE_LARGE.protein.k10.scaled5.sig.zip"
rm -rf "$TEMP_DIR/rust_output3.db"

# Rust version
echo "Rust CLI:"
RUST_START=$(date +%s.%N)
kmerseek-rust index --input "$TEST_FILE_LARGE" --output "$TEMP_DIR/rust_output3.db" --ksize 10 --encoding protein
RUST_END=$(date +%s.%N)
RUST_TIME=$(echo "$RUST_END - $RUST_START" | bc)
RUST_SIZE=$(du -h "$TEMP_DIR/rust_output3.db" | cut -f1)

# Python version
echo "Python CLI:"
PYTHON_START=$(date +%s.%N)
kmerseek index "$TEST_FILE_LARGE" --extract-kmers --moltype protein --ksize 10 --scaled 5 --force
PYTHON_END=$(date +%s.%N)
PYTHON_TIME=$(echo "$PYTHON_END - $PYTHON_START" | bc)
PYTHON_SIZE=$(du -h "$TEST_FILE_LARGE.protein.k10.scaled5.sig.zip" | cut -f1)

echo
echo "Results for large uncharacterized protein dataset:"
echo "  Rust:   ${RUST_TIME}s, ${RUST_SIZE}"
echo "  Python: ${PYTHON_TIME}s, ${PYTHON_SIZE}"
echo "  Speedup: $(echo "scale=2; $PYTHON_TIME / $RUST_TIME" | bc)x"
echo

echo "Test 4: Memory usage comparison"
echo "-------------------------------"

# Memory usage for larger file
echo "Rust memory usage:"
RUST_MEMORY=$(/usr/bin/time -l kmerseek-rust index --input "$TEST_FILE_LARGE" --output "$TEMP_DIR/rust_memory.db" --ksize 10 --encoding protein 2>&1 | grep "maximum resident set size" | awk '{print $1}')

echo "Python memory usage:"
PYTHON_MEMORY=$(/usr/bin/time -l kmerseek index "$TEST_FILE_LARGE" --extract-kmers --moltype protein --ksize 10 --scaled 5 --force 2>&1 | grep "maximum resident set size" | awk '{print $1}')

echo
echo "Memory usage (KB):"
echo "  Rust:   ${RUST_MEMORY}"
echo "  Python: ${PYTHON_MEMORY}"
echo "  Memory ratio: $(echo "scale=2; $PYTHON_MEMORY / $RUST_MEMORY" | bc)x"
echo

echo "=== Benchmark Complete ===" 