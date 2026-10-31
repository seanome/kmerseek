#!/bin/bash
#SBATCH --job-name=profile-chance-minimum
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH --output=profile_chance_minimum_%j.log
# Swiss-Prot search time, and how many of the unfiltered pairs are kept, at several
# --max-shared-kmer-pvalue cutoffs. Reuses the indexes and the 100 queries of
# scripts/profile_search_steps.sbatch, from the checkout named by INDEXES.
#   INDEXES=$SCRATCH/kmerseek-step-timing/profile_search_steps \
#       sbatch scripts/profile_chance_minimum.sbatch   # from the repository root
# Writes profile_chance_minimum/timing.tsv, steps.tsv and kept.tsv. The first search of each
# index warms the disk cache and is not used.
set -euo pipefail

module load python/3.12.1 rust/1.90.0 llvm/17.0.6 cmake/3.31.4 gcc/14.2.0
LLVM_ROOT=/share/software/user/open/llvm/17.0.6
CPATH=$(printf '%s' "${CPATH:-}" | tr ':' '\n' | { grep -v "^$LLVM_ROOT" || true; } | paste -sd: -)
export CPATH LIBCLANG_PATH=$LLVM_ROOT/lib CC=gcc CXX=g++ CARGO_HOME=$SCRATCH/.cargo
export RAYON_NUM_THREADS=$SLURM_CPUS_PER_TASK KMERSEEK_PROFILE=1

INDEXES=${INDEXES:?set INDEXES to the profile_search_steps output folder}
QUERIES=${QUERIES:-$SCRATCH/kmerseek-pr133/run_evalue_timing/queries.fasta}
OUT=$PWD/profile_chance_minimum
BIN=$PWD/target/release/kmerseek
mkdir -p "$OUT"

cargo build --release --bin kmerseek
git log --oneline -1 > "$OUT/commit.txt"

# alphabet, k, cutoff ("none" for no cutoff, "warmup" for the discarded first search)
RUNS=("protein20 10 warmup" "protein20 10 none" "protein20 10 1e-6"
      "hp_lehninger2 19 warmup" "hp_lehninger2 19 none" "hp_lehninger2 19 1e-2"
      "hp_lehninger2 19 1e-4" "hp_lehninger2 19 1e-6" "hp_lehninger2 19 1e-8")
printf 'alphabet\tk\tcutoff\twall_s\tcpu_s\tmax_rss_kb\tcsv_rows\n' > "$OUT/timing.tsv"
printf 'alphabet\tk\tcutoff\tstep\tseconds\tcalls\n' > "$OUT/steps.tsv"
for spec in "${RUNS[@]}"; do
    read -r alphabet k cutoff <<< "$spec"
    args=()
    case $cutoff in none|warmup) ;; *) args=(--max-shared-kmer-pvalue "$cutoff") ;; esac
    name=$alphabet.k$k.$cutoff
    /usr/bin/time -v -o "$OUT/time_$name.txt" "$BIN" search --query "$QUERIES" \
        --target "$INDEXES/swissprot.$alphabet.k$k.db" --output "$OUT/search_$name.csv" \
        --ksize "$k" --alphabet "$alphabet" ${args[@]+"${args[@]}"} > "$OUT/search_$name.log" 2>&1
    [ "$cutoff" = warmup ] && continue
    python3 - "$OUT/time_$name.txt" "$OUT/search_$name.csv" "$alphabet" "$k" "$cutoff" \
        >> "$OUT/timing.tsv" <<'PY'
import sys
t = dict(l.strip().rsplit(": ", 1) for l in open(sys.argv[1]) if ": " in l)
h, m, s = (["0"] + t["Elapsed (wall clock) time (h:mm:ss or m:ss)"].split(":"))[-3:]
wall = int(h) * 3600 + int(m) * 60 + float(s)
cpu = float(t["User time (seconds)"]) + float(t["System time (seconds)"])
rows = sum(1 for _ in open(sys.argv[2])) - 1
print("\t".join(sys.argv[3:6]), f"{wall:.2f}\t{cpu:.2f}\t{t['Maximum resident set size (kbytes)']}\t{rows}", sep="\t")
PY
    grep -P '^profile\t(?!step)' "$OUT/search_$name.log" | cut -f2- \
        | sed "s/^/$alphabet\t$k\t$cutoff\t/" >> "$OUT/steps.tsv" || true
done

# Share of the unfiltered pairs each cutoff keeps: all pairs, and the pairs whose best region
# score is in the top 1% and top 10% of the unfiltered run of the same index.
python3 - "$OUT" <<'PY' > "$OUT/kept.tsv"
import csv, sys
out = sys.argv[1]
def pairs(name):
    best = {}
    with open(f"{out}/search_{name}.csv") as f:
        for row in csv.DictReader(f):
            key = (row["query_md5"], row["target_md5"])
            best[key] = max(best.get(key, float("-inf")), float(row["region_poisson_score"]))
    return best
print("alphabet\tk\tcutoff\tpairs\tpairs_kept_fraction\ttop1pct_kept_fraction"
      "\ttop10pct_kept_fraction\ttop1pct_min_score\ttop10pct_min_score")
for alphabet, k, cutoffs in [("protein20", 10, ["1e-6"]),
                             ("hp_lehninger2", 19, ["1e-2", "1e-4", "1e-6", "1e-8"])]:
    full = pairs(f"{alphabet}.k{k}.none")
    scores = sorted(full.values())
    cut1, cut10 = scores[int(0.99 * len(scores))], scores[int(0.90 * len(scores))]
    top1 = [p for p, v in full.items() if v >= cut1]
    top10 = [p for p, v in full.items() if v >= cut10]
    for cutoff in ["none"] + cutoffs:
        kept = pairs(f"{alphabet}.k{k}.{cutoff}")
        frac = lambda keys: sum(p in kept for p in keys) / len(keys)
        print(alphabet, k, cutoff, len(kept), f"{len(kept) / len(full):.4f}", f"{frac(top1):.4f}",
              f"{frac(top10):.4f}", f"{cut1:.2f}", f"{cut10:.2f}", sep="\t")
PY
cat "$OUT/timing.tsv" "$OUT/kept.tsv"
