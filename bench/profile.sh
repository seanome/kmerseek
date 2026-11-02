#!/bin/bash
#SBATCH --job-name=kmerseek-profile
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --output=bench/profile_%j.log
# Where the time and memory of the bench/speed.sh kmerseek search go. Runs the same search
# twice on the index speed.sh built: once under `perf record` for a flame graph and the
# functions with the most self time, once under heaptrack for the allocation sites holding
# the most memory at peak.
#   INDEX=$SCRATCH/kmerseek-bench/speed_<job>/km.db sbatch bench/profile.sh
# Writes bench/baseline_flame.svg, bench/profile/perf_self.txt and
# bench/profile/heaptrack_peaks.txt.
set -euo pipefail

module load python/3.12.1 rust/1.90.0 llvm/17.0.6 cmake/3.31.4 gcc/14.2.0
LLVM_ROOT=/share/software/user/open/llvm/17.0.6
CPATH=$(printf '%s' "${CPATH:-}" | tr ':' '\n' | { grep -v "^$LLVM_ROOT" || true; } | paste -sd: -)
export CPATH LIBCLANG_PATH=$LLVM_ROOT/lib CC=gcc CXX=g++ CARGO_HOME=$SCRATCH/.cargo
export PATH=$CARGO_HOME/bin:$PATH RAYON_NUM_THREADS=8

: "${INDEX:?set INDEX to the km.db that bench/speed.sh built}"
TOOLS=/home/groups/ayeletv/opt/kmerseek-bench/env/bin
QUERIES=$PWD/bench/queries.fasta
BIN=$PWD/target/release/kmerseek
OUT=$SCRATCH/kmerseek-bench/profile_${SLURM_JOB_ID:-local}
REPORT=$PWD/bench/profile
mkdir -p "$OUT" "$REPORT"

cargo build --release --bin kmerseek
command -v inferno-flamegraph > /dev/null || cargo install inferno --locked
git rev-parse --short HEAD > "$REPORT/commit.txt"
search() {
    "$@" "$BIN" search --query "$QUERIES" --target "$INDEX" --ksize 24 \
        --alphabet hp_lehninger2 --extend-mismatch-penalty 2 --output "$OUT/km_hits.csv"
}

# CPU: sample user-space stacks 49 times a second per thread, unwound from DWARF line
# tables (the release profile keeps them).
search perf record -e cpu-clock -F 49 --call-graph dwarf -o "$OUT/perf.data" -- \
    > "$OUT/perf_search.log" 2>&1
perf script -i "$OUT/perf.data" 2> /dev/null | inferno-collapse-perf > "$OUT/stacks.folded"
inferno-flamegraph --title "kmerseek search, hp_lehninger2 k=24 scaled 10" \
    < "$OUT/stacks.folded" > "$PWD/bench/baseline_flame.svg"
perf report -i "$OUT/perf.data" --no-children --sort symbol --stdio --percent-limit 0.5 \
    2> /dev/null | grep -v '^$' > "$REPORT/perf_self.txt"

# Memory: every allocation with its call stack; the peaks list sums the bytes each
# allocation site held at the moment total heap use was highest.
search "$TOOLS/heaptrack" -o "$OUT/heaptrack" > "$OUT/heaptrack_search.log" 2>&1
"$TOOLS/heaptrack_print" -f "$(ls "$OUT"/heaptrack.*)" --print-peaks 1 --peak-limit 10 \
    --print-allocators 0 --print-temporary 0 --print-leaks 0 > "$REPORT/heaptrack_peaks.txt"
head -60 "$REPORT/perf_self.txt"
grep -A3 -m12 'peak' "$REPORT/heaptrack_peaks.txt" || head -80 "$REPORT/heaptrack_peaks.txt"
