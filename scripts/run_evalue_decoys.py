#!/usr/bin/env python3
"""Count chance runs on the decoy test and compare them with the run model.

Reproduces the numbers the README and `region_run_evalue_on_2mer_shuffled_decoys` quote
for runs seen against runs predicted: the 25 BCL-2-like test proteins against their 500
dipeptide-shuffled decoys, exact matching, split into a query against its own 20 shuffles
and a query against the 480 shuffles of other proteins.

A run is a maximal stretch of positions where query and decoy are in the same class, on
one diagonal (no gaps), at least k long. X (a residue outside the alphabet) never matches.
The model's count of runs at least L long in one pair is

    (1 - Pr(same)) (m - L + 1) (n - L + 1) Pr(same)^L

with Pr(same) from the two sequences' class compositions, m and n their lengths. This uses
each decoy's true length; kmerseek estimates the database's residue count from its k-mers,
so its E-values differ slightly (README, "E-values in the CSV").

    python scripts/run_evalue_decoys.py
"""

import gzip
from collections import Counter
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent.parent
QUERIES = HERE / "tests/testdata/fasta/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz"
DECOYS = HERE / "tests/testdata/fasta/bcl2_25_shuffled_2mer_20x_seed1.fasta.gz"

# hp_lehninger2 is sourmash's aa_to_hp; gbmr4 is GBMR4_CLUSTERS in src/rust/alphabets.rs.
ALPHABETS = {
    "hp_lehninger2": ["AFGILMPVWY", "NCSTDERHKQ"],
    "gbmr4": ["ADKERNTSQ", "YFLIVMCWH", "G", "P"],
    "protein20": list("ACDEFGHIKLMNPQRSTVWY"),
}
CASES = [("hp_lehninger2", 15), ("gbmr4", 10), ("protein20", 5)]


def read_fasta_gz(path):
    out, name, seq = [], None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(seq)))
                name, seq = line[1:], []
            elif line:
                seq.append(line)
    out.append((name, "".join(seq)))
    return out


def maximal_runs(q, db, k):
    """Length and database end of every maximal run of equal letters, at least k long,
    between q and db (several sequences joined by X), on any diagonal."""
    qa = np.frombuffer(q.encode(), dtype=np.uint8)
    da = np.frombuffer(db.encode(), dtype=np.uint8)
    x = ord("X")
    prev = np.zeros(len(da), dtype=np.int32)
    out = []
    for i in range(len(qa)):
        eq = (da == qa[i]) & (da != x) if qa[i] != x else np.zeros(len(da), bool)
        cur = np.zeros(len(da), dtype=np.int32)
        cur[1:] = np.where(eq[1:], prev[:-1] + 1, 0)
        cur[0] = int(eq[0])
        ended = np.nonzero((prev >= k) & ~np.r_[eq[1:], False])[0]
        out += [(j + 1, int(prev[j])) for j in ended]
        prev = cur
    out += [(j + 1, int(prev[j])) for j in np.nonzero(prev >= k)[0]]
    return out


def main():
    queries = read_fasta_gz(QUERIES)
    decoys = read_fasta_gz(DECOYS)
    # Each decoy is named "<source header>_shuffle<i>".
    source = [next(i for i, (qn, _) in enumerate(queries) if dn.startswith(qn + "_shuffle"))
              for dn, _ in decoys]
    for alphabet, k in CASES:
        table = {a: g[0].lower() for g in ALPHABETS[alphabet] for a in g}
        enc = lambda s: "".join(table.get(a, "X") for a in s.upper())
        letters = [g[0].lower() for g in ALPHABETS[alphabet]]
        qs = [enc(s) for _, s in queries]
        ts = [enc(s) for _, s in decoys]
        comp = lambda s: {c: s.count(c) / len(s) for c in letters}
        tc = [comp(t) for t in ts]
        db = "X".join(ts)
        starts = np.cumsum([0] + [len(t) + 1 for t in ts])
        lengths = [k, k + 2, k + 4]
        seen = {own: Counter() for own in (True, False)}
        pred = {own: Counter() for own in (True, False)}
        for qi, q in enumerate(qs):
            cq = comp(q)
            for ti, t in enumerate(ts):
                p = sum(cq[c] * tc[ti][c] for c in letters)
                for L in lengths:
                    pred[source[ti] == qi][L] += (
                        (1 - p) * max(len(q) - L + 1, 0) * max(len(t) - L + 1, 0) * p**L
                    )
            for d_end, run in maximal_runs(q, db, k):
                ti = int(np.searchsorted(starts, d_end - 1, side="right") - 1)
                for L in lengths:
                    seen[source[ti] == qi][L] += run >= L
        print(f"{alphabet} k={k}")
        for own, label in ((True, "a query vs its own 20 shuffles"),
                           (False, "a query vs the 480 shuffles of other proteins")):
            parts = [f"L>={L}: {seen[own][L]:_} seen, {pred[own][L]:_.0f} predicted" for L in lengths]
            print(f"  {label}: " + "; ".join(parts))


if __name__ == "__main__":
    main()
