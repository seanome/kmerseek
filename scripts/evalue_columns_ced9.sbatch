#!/bin/bash
#SBATCH --job-name=evalue-columns-ced9
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=evalue_columns_ced9_%j.log
# Makes the two CSVs behind README "Which one to use" and
# docs/images/evalue_columns_ced9_bcl2_human_k19.png: CED9 against the 300 human proteins
# of the BCL2 UniProt search, hp_lehninger2 k=19, once exact and once extended at
# --extend-mismatch-penalty 2. Run from the repository root:
#   sbatch scripts/evalue_columns_ced9.sbatch
# then draw the figure with scripts/plot_evalue_columns_ced9.py on ced9_k19_extended.csv.
set -euo pipefail

module load rust llvm/17.0.6 cmake gcc/14.2.0
# The llvm module puts libc++ headers on CPATH, which breaks the RocksDB build with g++.
CPATH=$(echo "${CPATH:-}" | tr ':' '\n' | { grep -v '/share/software/user/open/llvm/17.0.6' || true; } | paste -sd:)
export CPATH LIBCLANG_PATH=/share/software/user/open/llvm/17.0.6/lib CC=gcc CXX=g++

cargo build --release --no-default-features --bins
KMERSEEK=target/release/kmerseek
F=tests/testdata/fasta/uniprotkb_BCL2_AND_model_organism_9606_2025_02_06.fasta.gz
OUT=evalue_columns_ced9
mkdir -p "$OUT"

"$KMERSEEK" index --input "$F" --output "$OUT/idx.db" --ksize 19 --alphabet hp_lehninger2
for mode in exact extended; do
    extra=()
    [[ $mode == extended ]] && extra=(--extend-mismatch-penalty 2)
    "$KMERSEEK" search --query tests/testdata/fasta/ced9.fasta --target "$OUT/idx.db" \
        --output "$OUT/ced9_k19_$mode.csv" --ksize 19 --alphabet hp_lehninger2 \
        --min-shared-kmers 0 --max-pvalue 1 "${extra[@]}"
done
