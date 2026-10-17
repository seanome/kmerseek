#!/bin/bash
#SBATCH --job-name=ka-fit-curve
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=ka_fit_curve_%x_%j.log
# One fit curve under docs/images/ka_fit_curves/: index a FASTA at hp_thomas_dill2 k=12
# and fit r_database and K at C = 2, X = 8 on 200 calibration queries, seed 1, with one
# --ka-null. From the repository root, after `cargo build --release --bin kmerseek`:
#   sbatch --job-name=<name>_<null> scripts/ka_fit_curve.sbatch <fasta> <name> <null>
# <null> is the file-name spelling: database, shuffled, shuffled_dipeptide or reversed.
# Writes ka_fit_curves/<name>_<null>.csv and /usr/bin/time -v output next to it.
set -euo pipefail

FASTA=$1
NAME=$2
NULL=$3
OUT=$PWD/ka_fit_curves
KMERSEEK=$PWD/target/release/kmerseek
export RAYON_NUM_THREADS=$SLURM_CPUS_PER_TASK
mkdir -p "$OUT"
git log --oneline -1
"$KMERSEEK" --version

idx=$TMPDIR/${NAME}_${NULL}.db
/usr/bin/time -v -o "$OUT/${NAME}_${NULL}.time.txt" "$KMERSEEK" index --input "$FASTA" \
    --output "$idx" --ksize 12 --alphabet hp_thomas_dill2 \
    --extend-mismatch-penalty 2 --extend-xdrop 8 --ka-queries 200 --ka-seed 1 \
    --ka-null "${NULL//_/-}" --ka-survival-out "$OUT/${NAME}_${NULL}.csv"
command rm -rf "$idx"
grep -E "Elapsed|Maximum resident" "$OUT/${NAME}_${NULL}.time.txt"
