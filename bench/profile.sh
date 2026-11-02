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

# Peak resident memory split into heap and the like (RssAnon) and mapped files such as
# the index (RssFile), read from /proc every 0.1 s while the search runs.
"$BIN" search --query "$QUERIES" --target "$INDEX" --ksize 24 --alphabet hp_lehninger2 \
    --extend-mismatch-penalty 2 --output "$OUT/km_hits.csv" > "$OUT/rss_search.log" 2>&1 &
pid=$!
while kill -0 $pid 2> /dev/null; do
    { grep -E '^(VmRSS|RssAnon|RssFile|RssShmem):' /proc/$pid/status 2> /dev/null || true; } | tr '\n' ' '
    echo; sleep 0.1
done | awk '{ for (i = 1; i < NF; i += 3) if ($(i + 1) > max[$i]) max[$i] = $(i + 1) }
            END { for (k in max) printf "%s\tpeak %.2f GB\n", k, max[k] / 1024 / 1024 }' \
    | sort > "$REPORT/rss_split.txt"
wait $pid

# CPU: sample stacks 999 times a second per thread, unwound from DWARF line tables (the
# release profile keeps them). Kernel frames have no symbols here (perf_event_paranoid 2),
# so they show as [unknown]. Inlined frames are left out so time is charged to the
# function the compiler kept.
search perf record -e cpu-clock -F 999 --call-graph dwarf -o "$OUT/perf.data" -- \
    > "$OUT/perf_search.log" 2>&1
/usr/bin/time -v -o "$OUT/time_perf_script.txt" perf script --no-inline -i "$OUT/perf.data" \
    2> /dev/null > "$OUT/perf.script"
inferno-collapse-perf < "$OUT/perf.script" > "$OUT/stacks.folded"
inferno-flamegraph --title "kmerseek search, hp_lehninger2 k=24 scaled 10" \
    < "$OUT/stacks.folded" > "$PWD/bench/baseline_flame.svg"
python3 - "$OUT/stacks.folded" > "$REPORT/perf_self.txt" <<'PY2'
import sys
from collections import Counter
self_time, total = Counter(), 0
for line in open(sys.argv[1]):
    stack, n = line.rstrip().rsplit(" ", 1)
    self_time[stack.split(";")[-1]] += int(n)
    total += int(n)
print("share_of_samples\tfunction")
for name, n in self_time.most_common(30):
    print(f"{100 * n / total:.1f}%\t{name}")
PY2

# Memory: every allocation with its call stack; the peaks list sums the bytes each
# allocation site held at the moment total heap use was highest.
search "$TOOLS/heaptrack" -o "$OUT/heaptrack" > "$OUT/heaptrack_search.log" 2>&1
/usr/bin/time -v -o "$OUT/time_heaptrack_print.txt" "$TOOLS/heaptrack_print" \
    -f "$(ls "$OUT"/heaptrack.*)" --print-peaks 1 --peak-limit 10 --print-allocators 0 \
    --print-temporary 0 --print-leaks 0 > "$REPORT/heaptrack_peaks.txt"
grep -h 'Maximum resident' "$OUT"/time_*.txt
cat "$REPORT/rss_split.txt"
head -12 "$REPORT/perf_self.txt"
tail -4 "$REPORT/heaptrack_peaks.txt"
