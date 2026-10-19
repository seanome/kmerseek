#!/bin/bash
#SBATCH --job-name=run-evalue-timing
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=run_evalue_timing_%j.log
# How much time the run E-value adds to an exact search (PR 133 review): every candidate
# pair now counts the target's class composition (pair_pr_same) and walks each region
# again (measure_runs). Times kmerseek at main (before) and at this branch (after) on one
# node: 100 human proteins, picked with a fixed seed, searched exactly against the human
# reference proteome, hp_lehninger2 k=19, three runs per build, alternating. Each build
# gets its own index, since the index format differs (schema 4 and 5).
#   sbatch scripts/time_run_evalue_overhead.sbatch   # from the repository root
# Writes run_evalue_timing/timing.tsv: build, run, wall seconds, user+system CPU seconds,
# peak RSS in kB, and CSV rows written.
set -euo pipefail

module load python/3.12.1 rust/1.90.0 llvm/17.0.6 cmake/3.31.4 gcc/14.2.0
LLVM_ROOT=/share/software/user/open/llvm/17.0.6
CPATH=$(printf '%s' "${CPATH:-}" | tr ':' '\n' | { grep -v "^$LLVM_ROOT" || true; } | paste -sd: -)
export CPATH LIBCLANG_PATH=$LLVM_ROOT/lib CC=gcc CXX=g++ CARGO_HOME=$SCRATCH/.cargo
export RAYON_NUM_THREADS=$SLURM_CPUS_PER_TASK

PROTEOME=$SCRATCH/why-numbers/2026_03/inputs/UP000005640_9606.fasta.gz
OUT=$PWD/run_evalue_timing
BASE=${BASE:-origin/main}
mkdir -p "$OUT"

# Before: main, in its own clone. After: this checkout.
if [ ! -d "$OUT/main-src/.git" ]; then
    git clone --quiet https://github.com/seanome/kmerseek "$OUT/main-src"
fi
(cd "$OUT/main-src" && git fetch --quiet origin && git checkout --quiet "$BASE" && \
    git log --oneline -1 > "$OUT/main_commit.txt" && cargo build --release --bin kmerseek)
cargo build --release --bin kmerseek
git log --oneline -1 > "$OUT/pr_commit.txt"
declare -A BIN=([main]=$OUT/main-src/target/release/kmerseek [pr]=$PWD/target/release/kmerseek)

# 100 queries, the same for both builds.
python3 - "$PROTEOME" "$OUT/queries.fasta" <<'PY'
import gzip, random, sys
records, name, seq = [], None, []
for line in gzip.open(sys.argv[1], "rt"):
    if line.startswith(">"):
        if name: records.append((name, "".join(seq)))
        name, seq = line, []
    else:
        seq.append(line.strip())
records.append((name, "".join(seq)))
random.seed(20261017)
with open(sys.argv[2], "w") as out:
    for name, seq in random.sample(records, 100):
        out.write(f"{name}{seq}\n")
PY

for build in main pr; do
    "${BIN[$build]}" index --input "$PROTEOME" --output "$OUT/idx_$build.db" \
        --ksize 19 --alphabet hp_lehninger2 --ka-queries 0
done

printf 'build\trun\twall_s\tcpu_s\tmax_rss_kb\tcsv_rows\n' > "$OUT/timing.tsv"
for run in 1 2 3; do
    for build in main pr; do
        csv=$OUT/search_${build}_$run.csv
        /usr/bin/time -v -o "$OUT/time_${build}_$run.txt" "${BIN[$build]}" search \
            --query "$OUT/queries.fasta" --target "$OUT/idx_$build.db" --output "$csv" \
            --ksize 19 --alphabet hp_lehninger2 > "$OUT/search_${build}_$run.log" 2>&1
        python3 - "$OUT/time_${build}_$run.txt" "$csv" "$build" "$run" >> "$OUT/timing.tsv" <<'PY'
import sys
t = dict(l.strip().rsplit(": ", 1) for l in open(sys.argv[1]) if ": " in l)
h, m, s = (["0"] + t["Elapsed (wall clock) time (h:mm:ss or m:ss)"].split(":"))[-3:]
wall = int(h) * 3600 + int(m) * 60 + float(s)
cpu = float(t["User time (seconds)"]) + float(t["System time (seconds)"])
rows = sum(1 for _ in open(sys.argv[2])) - 1
print(f"{sys.argv[3]}\t{sys.argv[4]}\t{wall:.2f}\t{cpu:.2f}\t{t['Maximum resident set size (kbytes)']}\t{rows}")
PY
    done
done
cat "$OUT/timing.tsv"
