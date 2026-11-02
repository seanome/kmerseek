#!/bin/bash
#SBATCH --job-name=kmerseek-speed
#SBATCH --account=ayeletv
#SBATCH --partition=hns,normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=08:00:00
#SBATCH --output=bench/speed_%j.log
# Speed and memory of kmerseek, DIAMOND and MMseqs2 on one fixed benchmark: the 1,000
# Swiss-Prot queries in bench/queries.fasta against all of Swiss-Prot 2026_03, 8 threads
# each, one CPU node. The kmerseek setting is fixed too: hp_lehninger2, k=24, scaled 10,
# extension on (--extend-mismatch-penalty 2). See bench/README.md.
#   sbatch bench/speed.sh        # from the repository root
# Appends one row per tool and step to bench/ledger.tsv. Writes kmerseek's hit list,
# sorted, to $OUT/kmerseek_hits.tsv and compares it with bench/baseline_hits.tsv; with no
# baseline yet, it writes that file instead.
set -euo pipefail

module load python/3.12.1 rust/1.90.0 llvm/17.0.6 cmake/3.31.4 gcc/14.2.0
LLVM_ROOT=/share/software/user/open/llvm/17.0.6
CPATH=$(printf '%s' "${CPATH:-}" | tr ':' '\n' | { grep -v "^$LLVM_ROOT" || true; } | paste -sd: -)
export CPATH LIBCLANG_PATH=$LLVM_ROOT/lib CC=gcc CXX=g++ CARGO_HOME=$SCRATCH/.cargo

THREADS=8
export RAYON_NUM_THREADS=$THREADS
ALPHABET=hp_lehninger2
KSIZE=24
SCALED=10
PENALTY=2
TOOLS=/home/groups/ayeletv/opt/kmerseek-bench/env/bin
TARGET_GZ=${TARGET_GZ:-$SCRATCH/why-numbers/2026_03/inputs/uniprot_sprot.fasta.gz}
TARGET_MD5=bc9d398533e6df582b563c6c03093bd0
QUERIES=$PWD/bench/queries.fasta
QUERIES_MD5=a4f84f1f590c917616982eef85435b85
OUT=${OUT:-$SCRATCH/kmerseek-bench/speed_${SLURM_JOB_ID:-local}}
LEDGER=$PWD/bench/ledger.tsv
BASELINE=$PWD/bench/baseline_hits.tsv
BIN=$PWD/target/release/kmerseek

check_md5() {
    [ "$(md5sum < "$1" | cut -d' ' -f1)" = "$2" ] || { echo "md5 of $1 is not $2" >&2; exit 1; }
}
check_md5 "$TARGET_GZ" "$TARGET_MD5"
check_md5 "$QUERIES" "$QUERIES_MD5"

cargo build --release --bin kmerseek
COMMIT=$(git rev-parse --short HEAD)
DATE=$(date +%F)
CPU=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ *//; s/  */ /g')
NODE="$(hostname -s) $CPU"
mkdir -p "$OUT"
scontrol show node "$(hostname -s)" > "$OUT/node.txt" 2>/dev/null || true
{ "$TOOLS/diamond" version; "$TOOLS/mmseqs" version; "$BIN" --version; } > "$OUT/versions.txt"
[ -f "$LEDGER" ] || printf 'date\tgit_commit\ttool\tstep\twall_s\tpeak_rss_gb\tindex_bytes\tnode\tthreads\n' > "$LEDGER"
gzip -dc "$TARGET_GZ" > "$OUT/db.fasta"
cd "$OUT"

# Run a command under /usr/bin/time -v and append its ledger row. index_bytes is the
# on-disk size of the paths after the command, or "NA" when none are given.
measure() {
    local tool=$1 step=$2 sizes=$3; shift 3
    /usr/bin/time -v -o "time_${tool}_$step.txt" "$@" > "${tool}_$step.log" 2>&1
    local bytes=NA
    [ -n "$sizes" ] && bytes=$(du -scb $sizes | tail -1 | cut -f1)
    awk -F': ' -v pre="$DATE\t$COMMIT\t$tool\t$step" -v post="$bytes\t$NODE\t$THREADS" '
        /Elapsed \(wall clock\)/ { n = split($2, t, ":"); wall = 0
                                   for (i = 1; i <= n; i++) wall = wall * 60 + t[i] }
        /Maximum resident set size/ { rss = $2 / 1024 / 1024 }
        END { printf "%s\t%.2f\t%.2f\t%s\n", pre, wall, rss, post }
    ' "time_${tool}_$step.txt" | tee -a "$LEDGER"
}

measure kmerseek index "km.db" "$BIN" index --input db.fasta --output km.db \
    --ksize $KSIZE --alphabet $ALPHABET --scaled $SCALED --extend-mismatch-penalty $PENALTY
measure kmerseek search "" "$BIN" search --query "$QUERIES" --target km.db \
    --ksize $KSIZE --alphabet $ALPHABET --extend-mismatch-penalty $PENALTY \
    --output km_hits.csv

# The correctness gate: query, target, both regions and n_shared, one line per region,
# sorted bytewise so two runs compare with cmp. Python's csv module, because protein
# names contain commas and the CSV quotes them.
python3 - km_hits.csv > kmerseek_hits.unsorted.tsv <<'PY2'
import csv, sys
cols = ["query_name", "target_name", "region_start", "region_end", "target_start", "target_end",
        "region_n_shared_kmers"]
with open(sys.argv[1], newline="") as f:
    for row in csv.DictReader(f):
        sys.stdout.write("\t".join(row[c] for c in cols) + "\n")
PY2
{ printf 'query_name\ttarget_name\tregion_start\tregion_end\ttarget_start\ttarget_end\tn_shared\n'
  LC_ALL=C sort -S 16G --parallel=$THREADS kmerseek_hits.unsorted.tsv; } > kmerseek_hits.tsv
rm -f kmerseek_hits.unsorted.tsv
echo "kmerseek hit rows: $(($(wc -l < kmerseek_hits.tsv) - 1))"
if [ -f "$BASELINE" ]; then
    cmp kmerseek_hits.tsv "$BASELINE" && echo "HITS IDENTICAL to $BASELINE" \
        || echo "HITS DIFFER from $BASELINE"
else
    \cp kmerseek_hits.tsv "$BASELINE" && echo "wrote new baseline $BASELINE"
fi

measure diamond makedb "dmdb.dmnd" "$TOOLS/diamond" makedb --in db.fasta -d dmdb \
    -p $THREADS --quiet
measure diamond blastp "" "$TOOLS/diamond" blastp -q "$QUERIES" -d dmdb -o dm_hits.tsv \
    --sensitive -p $THREADS --quiet

measure mmseqs createdb "" "$TOOLS/mmseqs" createdb db.fasta mmdb -v 1
measure mmseqs createindex "mmdb*" "$TOOLS/mmseqs" createindex mmdb mmtmp \
    --threads $THREADS -v 1
measure mmseqs search "" "$TOOLS/mmseqs" easy-search "$QUERIES" mmdb mm_hits.tsv mmtmp2 \
    -s 5.7 --threads $THREADS -v 1

rm -f db.fasta
