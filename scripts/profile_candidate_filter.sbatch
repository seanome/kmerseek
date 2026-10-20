#!/bin/bash
#SBATCH --job-name=profile-candidate-filter
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH --output=profile_candidate_filter_%j.log
# Remeasures the Swiss-Prot search after candidates are filtered on their shared k-mer count
# before any target is read from disk, and tries --max-kmer-targets caps. Reuses the indexes
# and the 100 queries of scripts/profile_search_steps.sbatch (run that first).
#   sbatch scripts/profile_candidate_filter.sbatch   # from the repository root
# Writes profile_candidate_filter/timing.tsv, steps.tsv, and kept.tsv (how many of the uncapped
# results each cap keeps). Each search runs twice with the step timers on.
set -euo pipefail

module load python/3.12.1 rust/1.90.0 llvm/17.0.6 cmake/3.31.4 gcc/14.2.0
LLVM_ROOT=/share/software/user/open/llvm/17.0.6
CPATH=$(printf '%s' "${CPATH:-}" | tr ':' '\n' | { grep -v "^$LLVM_ROOT" || true; } | paste -sd: -)
export CPATH LIBCLANG_PATH=$LLVM_ROOT/lib CC=gcc CXX=g++ CARGO_HOME=$SCRATCH/.cargo
export RAYON_NUM_THREADS=$SLURM_CPUS_PER_TASK KMERSEEK_PROFILE=1

BEFORE=$PWD/profile_search_steps
QUERIES=${QUERIES:-$SCRATCH/kmerseek-pr133/run_evalue_timing/queries.fasta}
OUT=$PWD/profile_candidate_filter
BIN=$PWD/target/release/kmerseek
mkdir -p "$OUT"

cargo build --release --bin kmerseek
git log --oneline -1 > "$OUT/commit.txt"

# alphabet, k, cap ("none" for no cap)
RUNS=("protein20 10 none" "hp_lehninger2 19 none" "hp_lehninger2 19 10000"
      "hp_lehninger2 19 1000" "hp_lehninger2 19 100")
printf 'alphabet\tk\tcap\trun\twall_s\tcpu_s\tmax_rss_kb\tcsv_rows\n' > "$OUT/timing.tsv"
printf 'alphabet\tk\tcap\trun\tstep\tseconds\tcalls\n' > "$OUT/steps.tsv"
for spec in "${RUNS[@]}"; do
    read -r alphabet k cap <<< "$spec"
    cap_arg=(); [ "$cap" = none ] || cap_arg=(--max-kmer-targets "$cap")
    for run in 1 2; do
        name=$alphabet.k$k.cap$cap.run$run
        /usr/bin/time -v -o "$OUT/time_$name.txt" "$BIN" search --query "$QUERIES" \
            --target "$BEFORE/swissprot.$alphabet.k$k.db" --output "$OUT/search_$name.csv" \
            --ksize "$k" --alphabet "$alphabet" "${cap_arg[@]}" > "$OUT/search_$name.log" 2>&1
        python3 - "$OUT/time_$name.txt" "$OUT/search_$name.csv" \
            "$alphabet" "$k" "$cap" "$run" >> "$OUT/timing.tsv" <<'PY'
import sys
t = dict(l.strip().rsplit(": ", 1) for l in open(sys.argv[1]) if ": " in l)
h, m, s = (["0"] + t["Elapsed (wall clock) time (h:mm:ss or m:ss)"].split(":"))[-3:]
wall = int(h) * 3600 + int(m) * 60 + float(s)
cpu = float(t["User time (seconds)"]) + float(t["System time (seconds)"])
rows = sum(1 for _ in open(sys.argv[2])) - 1
print("\t".join(sys.argv[3:7]), f"{wall:.2f}\t{cpu:.2f}\t{t['Maximum resident set size (kbytes)']}\t{rows}", sep="\t")
PY
        grep -P '^profile\t(?!step)' "$OUT/search_$name.log" | cut -f2- \
            | sed "s/^/$alphabet\t$k\t$cap\t$run\t/" >> "$OUT/steps.tsv" || true
    done
done

# Without a cap the filter must not change a single row: compare with the run before it.
md5sum "$BEFORE"/search_protein20.k10.run2.csv "$OUT"/search_protein20.k10.capnone.run1.csv \
    "$BEFORE"/search_hp_lehninger2.k19.run2.csv "$BEFORE"/search_hp_lehninger2.k19.run4.csv \
    "$OUT"/search_hp_lehninger2.k19.capnone.run1.csv \
    > "$OUT/identical.md5"

# How many of the uncapped query-target pairs each cap keeps: all pairs, and the pairs whose
# best region score is in the top 1% and the top 10% of the uncapped run.
python3 - "$OUT" <<'PY' > "$OUT/kept.tsv"
import csv, sys
out = sys.argv[1]
def pairs(cap):
    best = {}
    with open(f"{out}/search_hp_lehninger2.k19.cap{cap}.run1.csv") as f:
        for row in csv.DictReader(f):
            key = (row["query_md5"], row["target_md5"])
            best[key] = max(best.get(key, float("-inf")), float(row["region_poisson_score"]))
    return best
full = pairs("none")
scores = sorted(full.values())
cut1, cut10 = scores[int(0.99 * len(scores))], scores[int(0.90 * len(scores))]
top1 = [k for k, v in full.items() if v >= cut1]
top10 = [k for k, v in full.items() if v >= cut10]
print("cap\tpairs\tpairs_kept_fraction\ttop1pct_kept_fraction\ttop10pct_kept_fraction"
      "\ttop1pct_min_score\ttop10pct_min_score")
for cap in ["none", "10000", "1000", "100"]:
    kept = pairs(cap)
    frac = lambda keys: sum(k in kept for k in keys) / len(keys)
    print(cap, len(kept), f"{len(kept) / len(full):.4f}", f"{frac(top1):.4f}",
          f"{frac(top10):.4f}", f"{cut1:.2f}", f"{cut10:.2f}", sep="\t")
PY
cat "$OUT/timing.tsv" "$OUT/identical.md5" "$OUT/kept.tsv"
