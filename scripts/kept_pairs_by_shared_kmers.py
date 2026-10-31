"""Pairs with the largest possible region score, and the top pairs ranked by their best
region's shared k-mer count, kept at each --max-shared-kmer-pvalue cutoff.

    python3 scripts/kept_pairs_by_shared_kmers.py profile_chance_minimum > kept_by_shared.tsv
"""

import csv
import sys

out = sys.argv[1]
def pairs(name):
    best = {}
    with open(f"{out}/search_{name}.csv") as f:
        for row in csv.DictReader(f):
            key = (row["query_md5"], row["target_md5"])
            score, shared = float(row["region_poisson_score"]), int(row["region_n_shared_kmers"])
            s0, n0 = best.get(key, (float("-inf"), 0))
            best[key] = (max(s0, score), max(n0, shared))
    return best
print("alphabet\tcutoff\tpairs\tmax_score_pairs\tmax_score_kept\ttop1pct_by_shared\ttop1pct_by_shared_kept\ttop1pct_shared_min\ttop10pct_by_shared_kept\ttop10pct_shared_min")
for alphabet, k, cutoffs in [("protein20", 10, ["1e-6"]), ("hp_lehninger2", 19, ["1e-2", "1e-4", "1e-6", "1e-8"])]:
    full = pairs(f"{alphabet}.k{k}.none")
    top_score = max(s for s, _ in full.values())
    maxed = [p for p, (s, _) in full.items() if s >= top_score]
    ranked = sorted(full, key=lambda p: -full[p][1])
    top1, top10 = ranked[: len(ranked) // 100], ranked[: len(ranked) // 10]
    for cutoff in ["none"] + cutoffs:
        kept = pairs(f"{alphabet}.k{k}.{cutoff}")
        frac = lambda keys: sum(p in kept for p in keys) / len(keys)
        print(alphabet, cutoff, len(kept), len(maxed), f"{frac(maxed):.4f}", len(top1), f"{frac(top1):.4f}",
              full[top1[-1]][1], f"{frac(top10):.4f}", full[top10[-1]][1], sep="\t")
