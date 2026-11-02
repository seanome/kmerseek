#!/bin/bash
#SBATCH --job-name=profile-posting-files
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --output=profile_posting_files_%j.log
# Searches Swiss-Prot (protein20 k=10, hp_lehninger2 k=19) before and after the inverted index
# moved from HashMaps loaded at startup to posting files read in place. "Before" is the binary
# and the schema 4 indexes of scripts/profile_search_steps.sbatch; "after" is this checkout,
# with the indexes rebuilt, since schema 5 cannot read schema 4. The two alternate, three runs
# each, on the same node, with the step timers on. The first run of each reads from a cold
# disk cache. Peak memory: hp_lehninger2 k=19 "before" reached 54 GB in job 48690134.
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
cargo build --release --bin kmerseek
AFTER_BIN=$PWD/target/release/kmerseek
{ git log --oneline -1; (cd "$BEFORE_DIR" && git log --oneline -1); } > "$OUT/commit.txt"
md5sum "$SWISSPROT" "$QUERIES" > "$OUT/inputs.md5"

SETTINGS=("protein20 10" "hp_lehninger2 19")
for setting in "${SETTINGS[@]}"; do
    read -r alphabet k <<< "$setting"
    db=$OUT/swissprot.$alphabet.k$k.db
    [ -f "$db/CURRENT" ] || /usr/bin/time -v -o "$OUT/time_index.$alphabet.k$k.txt" \
        "$AFTER_BIN" index --input "$SWISSPROT" --output "$db" --ksize "$k" \
        --alphabet "$alphabet" --ka-queries 0 > "$OUT/index.$alphabet.k$k.log" 2>&1
done

printf 'alphabet\tk\tversion\trun\twall_s\tcpu_s\tmax_rss_kb\tcsv_rows\n' > "$OUT/timing.tsv"
printf 'alphabet\tk\tversion\trun\tstep\tseconds\tcalls\n' > "$OUT/steps.tsv"
printf 'what\tbytes\n' > "$OUT/sizes.tsv"
for setting in "${SETTINGS[@]}"; do
    read -r alphabet k <<< "$setting"
    before_db=$BEFORE_DIR/profile_search_steps/swissprot.$alphabet.k$k.db
    after_db=$OUT/swissprot.$alphabet.k$k.db
    for run in 1 2 3; do
        for version in before after; do
            if [ "$version" = before ]; then bin=$OUT/kmerseek.before db=$before_db
            else bin=$AFTER_BIN db=$after_db; fi
            name=$alphabet.k$k.$version.run$run
            /usr/bin/time -v -o "$OUT/time_$name.txt" "$bin" search --query "$QUERIES" \
                --target "$db" --output "$OUT/search_$name.csv" \
                --ksize "$k" --alphabet "$alphabet" > "$OUT/search_$name.log" 2>&1
            python3 - "$OUT/time_$name.txt" "$OUT/search_$name.csv" "$alphabet" "$k" \
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
                | sed "s/^/$alphabet\t$k\t$version\t$run\t/" >> "$OUT/steps.tsv" || true
        done
    done
    # On-disk size of each index, and of the posting files alone.
    printf '%s k=%s before index\t%s\n' "$alphabet" "$k" "$(du -sb "$before_db" | cut -f1)" >> "$OUT/sizes.tsv"
    printf '%s k=%s after index\t%s\n' "$alphabet" "$k" "$(du -sb "$after_db" | cut -f1)" >> "$OUT/sizes.tsv"
    for f in "$after_db"/postings/*; do
        printf '%s k=%s after %s\t%s\n' "$alphabet" "$k" "$(basename "$f")" "$(stat -c %s "$f")" >> "$OUT/sizes.tsv"
    done
done

# The change must not alter a single result row.
(cd "$OUT" && md5sum search_*.csv) > "$OUT/identical.md5"
cat "$OUT/timing.tsv" "$OUT/sizes.tsv" "$OUT/identical.md5"
