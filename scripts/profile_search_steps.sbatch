#!/bin/bash
#SBATCH --job-name=profile-search-steps
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH --output=profile_search_steps_%j.log
# Where the time of one kmerseek search goes, step by step. Indexes all of Swiss-Prot
# (575,748 proteins) under two settings, protein20 k=10 and hp_lehninger2 k=19, then
# searches the same 100 human proteins against each index four times: twice with the step
# timers off (KMERSEEK_PROFILE unset) and twice with them on, alternating, so the timers'
# own cost shows up as a wall-time difference.
#   sbatch scripts/profile_search_steps.sbatch   # from the repository root
# Writes profile_search_steps/timing.tsv (wall, CPU, peak RSS per run) and
# profile_search_steps/steps.tsv (seconds per step, from the profiled runs).
set -euo pipefail

module load python/3.12.1 rust/1.90.0 llvm/17.0.6 cmake/3.31.4 gcc/14.2.0
LLVM_ROOT=/share/software/user/open/llvm/17.0.6
CPATH=$(printf '%s' "${CPATH:-}" | tr ':' '\n' | { grep -v "^$LLVM_ROOT" || true; } | paste -sd: -)
export CPATH LIBCLANG_PATH=$LLVM_ROOT/lib CC=gcc CXX=g++ CARGO_HOME=$SCRATCH/.cargo
export RAYON_NUM_THREADS=$SLURM_CPUS_PER_TASK

SWISSPROT=$SCRATCH/why-numbers/2026_03/inputs/uniprot_sprot.fasta.gz
# 100 human proteins picked with seed 20261017 by scripts/time_run_evalue_overhead.sbatch.
QUERIES=${QUERIES:-$SCRATCH/kmerseek-pr133/run_evalue_timing/queries.fasta}
OUT=$PWD/profile_search_steps
BIN=$PWD/target/release/kmerseek
mkdir -p "$OUT"

cargo build --release --bin kmerseek
git log --oneline -1 > "$OUT/commit.txt"
md5sum "$SWISSPROT" "$QUERIES" > "$OUT/inputs.md5"

SETTINGS=("protein20 10" "hp_lehninger2 19")
for setting in "${SETTINGS[@]}"; do
    read -r alphabet k <<< "$setting"
    idx=$OUT/swissprot.$alphabet.k$k.db
    [ -f "$idx/CURRENT" ] || /usr/bin/time -v -o "$OUT/time_index.$alphabet.k$k.txt" \
        "$BIN" index --input "$SWISSPROT" --output "$idx" --ksize "$k" --alphabet "$alphabet" \
        --ka-queries 0 > "$OUT/index.$alphabet.k$k.log" 2>&1
done

printf 'alphabet\tk\trun\tprofile\twall_s\tcpu_s\tmax_rss_kb\tcsv_rows\n' > "$OUT/timing.tsv"
printf 'alphabet\tk\trun\tstep\tseconds\tcalls\n' > "$OUT/steps.tsv"
for setting in "${SETTINGS[@]}"; do
    read -r alphabet k <<< "$setting"
    for run in 1 2 3 4; do
        name=$alphabet.k$k.run$run
        if [ $((run % 2)) -eq 0 ]; then profile=on; export KMERSEEK_PROFILE=1
        else profile=off; unset KMERSEEK_PROFILE; fi
        /usr/bin/time -v -o "$OUT/time_$name.txt" "$BIN" search --query "$QUERIES" \
            --target "$OUT/swissprot.$alphabet.k$k.db" --output "$OUT/search_$name.csv" \
            --ksize "$k" --alphabet "$alphabet" > "$OUT/search_$name.log" 2>&1
        python3 - "$OUT/time_$name.txt" "$OUT/search_$name.csv" \
            "$alphabet" "$k" "$run" "$profile" >> "$OUT/timing.tsv" <<'PY'
import sys
t = dict(l.strip().rsplit(": ", 1) for l in open(sys.argv[1]) if ": " in l)
h, m, s = (["0"] + t["Elapsed (wall clock) time (h:mm:ss or m:ss)"].split(":"))[-3:]
wall = int(h) * 3600 + int(m) * 60 + float(s)
cpu = float(t["User time (seconds)"]) + float(t["System time (seconds)"])
rows = sum(1 for _ in open(sys.argv[2])) - 1
print("\t".join(sys.argv[3:7]), f"{wall:.2f}\t{cpu:.2f}\t{t['Maximum resident set size (kbytes)']}\t{rows}", sep="\t")
PY
        grep -P '^profile\t(?!step)' "$OUT/search_$name.log" | cut -f2- \
            | sed "s/^/$alphabet\t$k\t$run\t/" >> "$OUT/steps.tsv" || true
    done
done
unset KMERSEEK_PROFILE
cat "$OUT/timing.tsv" "$OUT/steps.tsv"
