#!/usr/bin/env python3
"""The two checks of `region_run_evalue` in src/rust/tests/test_run_evalue.rs, drawn.

    cargo test --no-default-features --lib write_run_evalue_figure_data -- --ignored
    python scripts/plot_run_evalue_checks.py

(a) Random sequences. Ten 1,000-residue queries against forty 1,000-residue targets, each
residue drawn on its own from the Swiss-Prot composition. For each minimum run length L0,
the runs an exact search reports are compared with the summed run E-value: (counted -
predicted) / sqrt(predicted), the distance in Poisson standard deviations. The formula
counts m x n start positions, but a run of L0 cannot start in the last L0 - 1 positions of
either sequence, so it predicts (1 - (L0 - 1) / 1000)^2 too many; the line shows the
distance that alone would give. The test asserts the filled points.

(b) Decoys. The 25 BCL-2-like test proteins against 500 decoys (each protein shuffled 20
times keeping its dipeptide counts). Query-target pairs whose best region has a run
E-value at or below each cutoff. A calibrated E-value would give, for each query, x chance
hits at E <= x over the whole database, spread over 500 targets: 12,500 (1 - e^(-x/500))
pairs, about 25 x while x is small.
"""

import argparse
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pubfig as pf  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

DATA = "docs/images/run_evalue_checks"
SEQUENCE_LENGTH = 1000
N_PAIRS = 25 * 500
N_TARGETS = 500
# L0 the test asserts, from region_run_evalue_matches_runs_between_random_sequences.
ASSERTED_L0 = {"hp_lehninger2": {18, 19, 20}, "gbmr4": {15, 16, 17}}
# Cutoffs the decoy test pins, from region_run_evalue_on_2mer_shuffled_decoys.
PINNED_CUTOFFS = (1, 10)

# One colour per alphabet-ksize pair, the same in both panels.
SERIES = {
    ("hp_lehninger2", 12): (pf.OKABE_ITO["blue"], "hp_lehninger2, k = 12 (h AFGILMPVWY, p CDEHKNQRST)"),
    ("gbmr4", 10): (pf.OKABE_ITO["vermillion"], "gbmr4, k = 10 (ADKERNTSQ, YFLIVMCWH, G, P)"),
    ("hp_lehninger2", 15): (pf.OKABE_ITO["bluish_green"], "hp_lehninger2, k = 15"),
}


def read(name):
    with open(os.path.join(DATA, name)) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def plot_random_sequences(ax):
    rows = read("random_sequences_runs.tsv")
    ax.axhspan(-3, 3, color=pf.GREY, alpha=0.2, lw=0, label="within 3 standard deviations")
    ax.axhline(0, color="black", lw=0.5)
    for (alphabet, ksize), (colour, label) in SERIES.items():
        mine = [r for r in rows if (r["alphabet"], int(r["ksize"])) == (alphabet, ksize)]
        if not mine:
            continue
        l0 = [int(r["l0"]) for r in mine]
        predicted = [float(r["predicted"]) for r in mine]
        z = [(int(r["counted"]) - p) / math.sqrt(p) for r, p in zip(mine, predicted)]
        edge = [p * ((1 - (x - 1) / SEQUENCE_LENGTH) ** 2 - 1) / math.sqrt(p) for x, p in zip(l0, predicted)]
        ax.plot(l0, edge, "--", color=colour, lw=0.8)
        asserted = [x in ASSERTED_L0[alphabet] for x in l0]
        ax.plot([x for x, a in zip(l0, asserted) if a], [v for v, a in zip(z, asserted) if a],
                "o", color=colour, ms=4)
        ax.plot([x for x, a in zip(l0, asserted) if not a], [v for v, a in zip(z, asserted) if not a],
                "o", mfc="white", color=colour, ms=4)
    ax.set_xlabel("minimum run length L0 (encoded positions)")
    ax.set_ylabel("(counted - predicted) / sqrt(predicted)\n(Poisson standard deviations)")
    ax.set_xticks(range(10, 23, 2))


def plot_decoys(ax):
    rows = read("decoys_pairs_at_or_below.tsv")
    xs = [10 ** (i / 20) for i in range(-10, 61)]
    ax.plot(xs, [N_PAIRS * (1 - math.exp(-x / N_TARGETS)) for x in xs], color="black", lw=1,
            label="calibrated E-value: 12,500 (1 - e^(-x/500))")
    for (alphabet, ksize), (colour, _) in SERIES.items():
        mine = [r for r in rows if (r["alphabet"], int(r["ksize"])) == (alphabet, ksize) and int(r["pairs"]) > 0]
        if mine:
            # k = 12 is drawn first, wide, so k = 15 on top of it leaves it visible.
            wide = ksize == 12
            ax.plot([float(r["evalue_cutoff"]) for r in mine], [int(r["pairs"]) for r in mine],
                    "o-", color=colour, ms=6 if wide else 2.5, lw=3 if wide else 1, zorder=2 if wide else 3)
    for x in PINNED_CUTOFFS:
        ax.axvline(x, color=pf.GREY, lw=0.5, ls=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([0.3, 1, 10, 100, 1000], ["0.3", "1", "10", "100", "1000"])
    ax.set_yticks([1, 10, 100, 1000, 10000], ["1", "10", "100", "1000", "10,000"])
    ax.set_xlim(0.3, 1000)
    ax.set_xlabel("region_run_evalue cutoff x\n(dotted: the cutoffs the test pins, 1 and 10)")
    ax.set_ylabel("query-target pairs with best region at E ≤ x\n(of 12,500)")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", default="docs/images/run_evalue_checks")
    args = parser.parse_args()

    pf.use_style()
    fig, (a, b) = pf.figure(pf.TWO_COLUMN_MM, 75, ncols=2)
    plot_random_sequences(a)
    plot_decoys(b)
    handles = [
        Line2D([], [], color=c, marker="o", ls="-", ms=4, label=label)
        for (_, k), (c, label) in SERIES.items()
    ]
    handles += [
        Line2D([], [], color="black", marker="o", mfc="white", ls="", ms=4,
                      label="(a) open: an L0 the test does not assert"),
        Line2D([], [], color="black", ls="--", lw=0.8,
                      label="(a) dashed: shortfall the sequence ends alone give"),
        Rectangle((0, 0), 1, 1, color=pf.GREY, alpha=0.2, label="(a) within 3 standard deviations"),
        Line2D([], [], color="black", lw=1, label="(b) calibrated E-value: 12,500 (1 - e^(-x/500))"),
    ]
    fig.legend(handles=handles, loc="outside upper center", ncol=2, frameon=False)
    pf.panel_label(a, "a")
    pf.panel_label(b, "b")
    pf.save(fig, args.output, formats=("png",))


if __name__ == "__main__":
    main()
