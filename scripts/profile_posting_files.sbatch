#!/bin/bash
#SBATCH --job-name=profile-posting-files
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=profile_posting_files_%j.log
# Searches Swiss-Prot (protein20, k=10) before and after the inverted index moved from
# HashMaps loaded at startup to k-mer-to-protein index files read in place. "Before" is the binary and the
# schema 4 index of scripts/profile_search_steps.sbatch; "after" is this checkout, with the
# index rebuilt, since schema 5 cannot read schema 4. The two alternate, three runs each, on
# the same node, with the step timers on. The first run of each reads from a cold disk cache.
#   sbatch scripts/profile_posting_files.sbatch   # from the repository root
# Writes profile_posting_files/timing.tsv, steps.tsv, sizes.tsv and identical.md5.
set -euo pipefail

module load python/3.12.1 rust/1.90.0 llvm/17.0.6 cmake/3.31.4 gcc/14.2.0
LLVM_ROOT=/share/software/user/open/llvm/17.0.6
CPATH=$(printf '%s' "${CPATH:-}" | tr ':' '\n' | { grep -v "^$LLVM_ROOT" || true; } | paste -sd: -)
export CPATH LIBCLANG_PATH=$LLVM_ROOT/lib CC=gcc CXX=g++ CARGO_HOME=$SCRATCH/.cargo
export RAYON_NUM_THREADS=$SLURM_CPUS_PER_TASK KMERSEEK_PROFILE=1

SWISSPROT=$SCRATCH/why-numbers/2026_03/inputs/uniprot_sprot.fasta.gz
QUERIES=${QUERIES:-$SCRATCH/kmerseek-pr133/run_evalue_timing/queries.fasta}
BEFORE_DIR=${BEFORE_DIR:-$SCRATCH/kmerseek-step-timing}
OUT=$PWD/profile_posting_files
mkdir -p "$OUT"

# Copied, so a later build in the other clone cannot change what "before" means.
/bin/cp -f "$BEFORE_DIR/target/release/kmerseek" "$OUT/kmerseek.before"
BEFORE_DB=$BEFORE_DIR/profile_search_steps/swissprot.protein20.k10.db
cargo build --release --bin kmerseek
AFTER_BIN=$PWD/target/release/kmerseek
AFTER_DB=$OUT/swissprot.protein20.k10.db
{ git log --oneline -1; (cd "$BEFORE_DIR" && git log --oneline -1); } > "$OUT/commit.txt"
md5sum "$SWISSPROT" "$QUERIES" > "$OUT/inputs.md5"

[ -f "$AFTER_DB/CURRENT" ] || /usr/bin/time -v -o "$OUT/time_index.txt" \
    "$AFTER_BIN" index --input "$SWISSPROT" --output "$AFTER_DB" --ksize 10 \
    --alphabet protein20 --ka-queries 0 > "$OUT/index.log" 2>&1

printf 'version\trun\twall_s\tcpu_s\tmax_rss_kb\tcsv_rows\n' > "$OUT/timing.tsv"
printf 'version\trun\tstep\tseconds\tcalls\n' > "$OUT/steps.tsv"
for run in 1 2 3; do
    for version in before after; do
        if [ "$version" = before ]; then bin=$OUT/kmerseek.before db=$BEFORE_DB
        else bin=$AFTER_BIN db=$AFTER_DB; fi
        name=$version.run$run
        /usr/bin/time -v -o "$OUT/time_$name.txt" "$bin" search --query "$QUERIES" \
            --target "$db" --output "$OUT/search_$name.csv" \
            --ksize 10 --alphabet protein20 > "$OUT/search_$name.log" 2>&1
        python3 - "$OUT/time_$name.txt" "$OUT/search_$name.csv" "$version" "$run" \
            >> "$OUT/timing.tsv" <<'PY'
import sys
t = dict(l.strip().rsplit(": ", 1) for l in open(sys.argv[1]) if ": " in l)
h, m, s = (["0"] + t["Elapsed (wall clock) time (h:mm:ss or m:ss)"].split(":"))[-3:]
wall = int(h) * 3600 + int(m) * 60 + float(s)
cpu = float(t["User time (seconds)"]) + float(t["System time (seconds)"])
rows = sum(1 for _ in open(sys.argv[2])) - 1
print(sys.argv[3], sys.argv[4], f"{wall:.2f}\t{cpu:.2f}\t{t['Maximum resident set size (kbytes)']}\t{rows}", sep="\t")
PY
        grep -P '^profile\t(?!step)' "$OUT/search_$name.log" | cut -f2- \
            | sed "s/^/$version\t$run\t/" >> "$OUT/steps.tsv" || true
    done
done

# On-disk size of each index, and of the k-mer-to-protein index files alone.
printf 'what\tbytes\n' > "$OUT/sizes.tsv"
printf 'before index\t%s\n' "$(du -sb "$BEFORE_DB" | cut -f1)" >> "$OUT/sizes.tsv"
printf 'after index\t%s\n' "$(du -sb "$AFTER_DB" | cut -f1)" >> "$OUT/sizes.tsv"
for f in "$AFTER_DB"/postings/*; do
    printf 'after %s\t%s\n' "$(basename "$f")" "$(stat -c %s "$f")" >> "$OUT/sizes.tsv"
done

# The change must not alter a single result row.
(cd "$OUT" && md5sum search_*.csv) > "$OUT/identical.md5"
cat "$OUT/timing.tsv" "$OUT/sizes.tsv" "$OUT/identical.md5"
