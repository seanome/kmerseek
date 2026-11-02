#!/bin/bash
#SBATCH --job-name=profile-query-frequencies
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=02:30:00
#SBATCH --output=profile_query_frequencies_%j.log
# Searches Swiss-Prot with hp_lehninger2 k=19 using three binaries, alternating, three runs
# each, on one node, with the step timers on:
#   before    the PR 158 binary on its schema 4 index (as in profile_posting_files.sbatch)
#   previous  posting files read in place, a frequency lookup in them for every shared k-mer
#             of every pair (the binary profile_posting_files.sbatch built, commit 3cddf3b)
#   after     this checkout: each query k-mer's frequency looked up once per query
# "previous" and "after" read the same schema 5 index that profile_posting_files.sbatch
# built. Run that script first. Peak memory: 55 GB in job 50498949; each run takes 3-8 min.
#   sbatch scripts/profile_query_frequencies.sbatch   # from the repository root
# Writes profile_query_frequencies/timing.tsv, steps.tsv and identical.md5.
set -euo pipefail

module load python/3.12.1 rust/1.90.0 llvm/17.0.6 cmake/3.31.4 gcc/14.2.0
LLVM_ROOT=/share/software/user/open/llvm/17.0.6
CPATH=$(printf '%s' "${CPATH:-}" | tr ':' '\n' | { grep -v "^$LLVM_ROOT" || true; } | paste -sd: -)
export CPATH LIBCLANG_PATH=$LLVM_ROOT/lib CC=gcc CXX=g++ CARGO_HOME=$SCRATCH/.cargo
export RAYON_NUM_THREADS=$SLURM_CPUS_PER_TASK KMERSEEK_PROFILE=1

SWISSPROT=$SCRATCH/why-numbers/2026_03/inputs/uniprot_sprot.fasta.gz
QUERIES=${QUERIES:-$SCRATCH/kmerseek-pr133/run_evalue_timing/queries.fasta}
BEFORE_DIR=${BEFORE_DIR:-$SCRATCH/kmerseek-step-timing}
ALPHABET=hp_lehninger2 K=19
POSTINGS_OUT=$PWD/profile_posting_files
OUT=$PWD/profile_query_frequencies
mkdir -p "$OUT"

# Copied before the build, so the build cannot change what "before" and "previous" mean.
[ -f "$OUT/kmerseek.previous" ] || /bin/cp -f "$PWD/target/release/kmerseek" "$OUT/kmerseek.previous"
/bin/cp -f "$BEFORE_DIR/target/release/kmerseek" "$OUT/kmerseek.before"
cargo build --release --bin kmerseek
/bin/cp -f "$PWD/target/release/kmerseek" "$OUT/kmerseek.after"
{ git log --oneline -1; (cd "$BEFORE_DIR" && git log --oneline -1); } > "$OUT/commit.txt"
md5sum "$OUT"/kmerseek.* "$SWISSPROT" "$QUERIES" > "$OUT/inputs.md5"

printf 'alphabet\tk\tversion\trun\twall_s\tcpu_s\tmax_rss_kb\tcsv_rows\n' > "$OUT/timing.tsv"
printf 'alphabet\tk\tversion\trun\tstep\tseconds\tcalls\n' > "$OUT/steps.tsv"
for run in 1 2 3; do
    for version in before previous after; do
        if [ "$version" = before ]; then db=$BEFORE_DIR/profile_search_steps/swissprot.$ALPHABET.k$K.db
        else db=$POSTINGS_OUT/swissprot.$ALPHABET.k$K.db; fi
        name=$ALPHABET.k$K.$version.run$run
        /usr/bin/time -v -o "$OUT/time_$name.txt" "$OUT/kmerseek.$version" search \
            --query "$QUERIES" --target "$db" --output "$OUT/search_$name.csv" \
            --ksize "$K" --alphabet "$ALPHABET" > "$OUT/search_$name.log" 2>&1
        python3 - "$OUT/time_$name.txt" "$OUT/search_$name.csv" "$ALPHABET" "$K" \
            "$version" "$run" >> "$OUT/timing.tsv" <<'PY'
import sys
t = dict(l.strip().rsplit(": ", 1) for l in open(sys.argv[1]) if ": " in l)
h, m, s = (["0"] + t["Elapsed (wall clock) time (h:mm:ss or m:ss)"].split(":"))[-3:]
wall = int(h) * 3600 + int(m) * 60 + float(s)
cpu = float(t["User time (seconds)"]) + float(t["System time (seconds)"])
rows = sum(1 for _ in open(sys.argv[2])) - 1
print(*sys.argv[3:7], f"{wall:.2f}\t{cpu:.2f}\t{t['Maximum resident set size (kbytes)']}\t{rows}", sep="\t")
PY
        grep -P '^profile\t(?!step)' "$OUT/search_$name.log" | cut -f2- \
            | sed "s/^/$ALPHABET\t$K\t$version\t$run\t/" >> "$OUT/steps.tsv" || true
    done
done

# The change must not alter a single result row.
(cd "$OUT" && md5sum search_*.csv) > "$OUT/identical.md5"
cat "$OUT/timing.tsv" "$OUT/identical.md5"
