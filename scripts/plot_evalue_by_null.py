#!/usr/bin/env python3
"""E-value against region score for fits read from different calibration queries, on one
database, so the effect of `--ka-null` on a reported E-value is visible.

    python scripts/plot_evalue_by_null.py -o docs/images/evalue_by_null_swissprot.png

E = K m n e^(-r_database x) for a 200-residue query (m) against all of Swiss-Prot
(n = 210 M residues), with r_database and K from the 15,000-sequence Swiss-Prot sample
fitted three ways (the curves under docs/images/ka_fit_curves).
"""

import argparse
import csv
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CURVES = "docs/images/ka_fit_curves/swissprot_sample15000_{}.csv"
QUERY_RESIDUES = 200
DATABASE_RESIDUES = 210e6
CUTOFF = 0.01
FITS = [
    ("database", "database sequences as they are (--ka-null database, the default)", "#1f77b4", "-"),
    ("shuffled_dipeptide", "shuffled keeping dipeptides (--ka-null shuffled-dipeptide)", "#7f7f7f", "--"),
    ("shuffled", "shuffled (--ka-null shuffled)", "#7f7f7f", ":"),
]


def read_fit(null):
    row = next(csv.DictReader(open(CURVES.format(null))))
    return float(row["r_database"]), float(row["k"])


def evalue(r, k, x):
    return k * QUERY_RESIDUES * DATABASE_RESIDUES * math.exp(-r * x)


def score_at(r, k, e):
    return math.log(k * QUERY_RESIDUES * DATABASE_RESIDUES / e) / r


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    xs = [10 + i * 0.1 for i in range(221)]
    fig, ax = plt.subplots(figsize=(10.5, 5.6), dpi=150)
    for null, label, colour, style in FITS:
        r, k = read_fit(null)
        x_cut = score_at(r, k, CUTOFF)
        ax.plot(xs, [evalue(r, k, x) for x in xs], style, color=colour, lw=2,
                label=f"{label}: r_database {r:.3f}, K {k:.4f}; E = {CUTOFF:g} at x = {x_cut:.1f}")
        ax.plot([x_cut], [CUTOFF], "o", color=colour, ms=6)
    ax.axhline(CUTOFF, color="#111111", ls="-", lw=0.6, label=f"E = {CUTOFF:g}, a common cutoff; the dots mark where each line crosses it")
    ax.set_yscale("log")
    ax.set_xlabel("region score x = λ_pair · S (nats)")
    ax.set_ylabel("E-value: regions this good expected by chance\n(200-residue query, all of Swiss-Prot)")
    ax.set_ylim(1e-5, 1e6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), frameon=False, fontsize=9)
    fig.suptitle("A fit on shuffled sequences gives the same region a smaller E-value\n(2 to 4 times smaller when dipeptides are kept, x = 15 to 25)",
                 x=0.02, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(args.output, bbox_inches="tight")


if __name__ == "__main__":
    main()
