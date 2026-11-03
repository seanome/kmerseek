#!/usr/bin/env python3
"""Draw what the protein20 calibration fix changes, on the 25 BCL2-family proteins of the
test fixture at protein20 k=5, mismatch penalty 0.14, give-up margin 0.57.

Panel a: the database's match probability (the chance two random positions hold the same
amino acid). Before the fix kmerseek counted it on the encoded sequence, which protein20
does not store, and got 0. After the fix it is the sum of squared residue frequencies,
counted here straight from the FASTA file as a check.

Panel b: the calibration regions' scores, lambda S in nats. Before the fix no protein20
region got a score, so all 64 sat at 0 and the fit was refused.

The kmerseek values are the ones
`cargo test --no-default-features --lib protein20` asserts
(`test_calibration_match_probability_protein20_counts_residues`,
`test_calibration_bins_protein20_spread_past_zero`).

    python scripts/plot_protein20_calibration.py \
        --output docs/images/protein20_calibration_bcl2_first25_k5
"""

import argparse
import gzip
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import pubfig as pf  # noqa: E402

FIXTURE = Path(__file__).parent.parent / (
    "tests/testdata/fasta/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz"
)
BEFORE = pf.OKABE_ITO["vermillion"]
AFTER = pf.OKABE_ITO["blue"]
BIN_WIDTH = 0.5  # nats, `evalue::BIN_WIDTH`

# kmerseek's match probability before and after the fix (from the tests).
KMERSEEK_BEFORE = 0.0
KMERSEEK_AFTER = 0.060_755_021_065_774_63
# Calibration bins (bin index, regions) from the tests: 25 shuffled queries, 64 regions.
BINS_BEFORE = {0: 64}
BINS_AFTER = {6: 2, 7: 1, 8: 7, 9: 7, 10: 16, 11: 12, 12: 12, 13: 5, 15: 1, 16: 1}


def residue_match_probability(path):
    """Sum of squared residue frequencies over every sequence in a gzipped FASTA."""
    counts = Counter()
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith(">"):
                counts.update(line.strip())
    total = sum(counts.values())
    return sum((n / total) ** 2 for n in counts.values()), total


def draw_match_probability(ax, by_hand):
    labels = ["kmerseek\nbefore fix", "kmerseek\nafter fix", "counted from\nthe FASTA"]
    values = [KMERSEEK_BEFORE, KMERSEEK_AFTER, by_hand]
    ax.bar(range(3), values, color=[BEFORE, AFTER, pf.GREY], width=0.6)
    for x, v in enumerate(values):
        ax.text(x, v + 0.002, f"{v:.3f}", ha="center", va="bottom")
    ax.set_xticks(range(3), labels)
    ax.set_ylim(0, 0.075)
    ax.set_ylabel("Match probability")


def draw_bins(ax):
    for bins, colour, label in [
        (BINS_BEFORE, BEFORE, "before fix"),
        (BINS_AFTER, AFTER, "after fix"),
    ]:
        xs = [b * BIN_WIDTH for b in bins]
        ax.bar(xs, list(bins.values()), width=BIN_WIDTH, align="edge", color=colour, label=label)
    ax.set_xlabel("Region score, lambda S (nats)")
    ax.set_ylabel("Calibration regions")
    ax.set_xlim(0, 9)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--output", required=True, help="path stem, no extension")
    args = parser.parse_args()
    by_hand, total = residue_match_probability(FIXTURE)
    assert total == 9288, total
    assert abs(by_hand - KMERSEEK_AFTER) < 1e-12, by_hand
    pf.use_style()
    fig, (ax_a, ax_b) = pf.figure(pf.ONE_COLUMN_MM, 55, ncols=2)
    draw_match_probability(ax_a, by_hand)
    draw_bins(ax_b)
    pf.shared_legend(fig)
    for ax, letter in zip((ax_a, ax_b), "ab"):
        pf.panel_label(ax, letter)
    pf.save(fig, args.output)


if __name__ == "__main__":
    main()
