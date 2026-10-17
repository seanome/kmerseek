#!/usr/bin/env python3
"""Draw one query against two database entries whose 5-mer sets are equal but whose
sequences differ: a poly-A run of 5 ("short") and of 8 ("long"), at protein20 k=5.

A poly-A run of any length from 5 up adds only one 5-mer, AAAAA, so both entries hold
the same 19 5-mers. The index used to key an entry by the sum of its k-mer hashes, so
"long" got the same key as "short" and was dropped. Keyed by its sequence, each entry is
its own target, and each region is cut from that target's own residues.

The query has a run of 7, so its key matches neither entry. The regions are the ones
`cargo test --no-default-features --lib same_kmer_set_different_sequence_is_its_own_target`
asserts, as (query start, query end, target start, target end), 0-based, end exclusive.

    python scripts/plot_same_kmer_set_targets.py \
        --output docs/images/same_kmer_set_polyA5_polyA8_protein20_k5
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import matplotlib.pyplot as plt  # noqa: E402
import pubfig as pf  # noqa: E402

QUERY = "MKTWLRDEHPAAAAAAACYFNQSGV"
TARGETS = {
    "short": "MKTWLRDEHPAAAAACYFNQSGV",
    "long": "MKTWLRDEHPAAAAAAAACYFNQSGV",
}
REGIONS = {
    "short": [(0, 15, 0, 15), (11, 16, 10, 15), (12, 25, 10, 23)],
    "long": [
        (0, 17, 0, 17),
        (10, 25, 11, 26),
        (10, 16, 12, 18),
        (10, 15, 13, 18),
        (11, 17, 10, 16),
        (12, 17, 10, 15),
    ],
}
FLANK = pf.OKABE_ITO["blue"]
POLY_A = pf.OKABE_ITO["vermillion"]
RUN_FILL = "#eeeeee"  # the poly-A run, behind every other mark
K = 5
CONCLUSION = (
    "Keyed by its sequence, each entry is its own target: search finds 3 regions in the "
    "short entry and 6 in the long one.\nBefore PR 152, the long entry was not stored "
    "and this search found no hits."
)


def kmer_set(seq):
    """The distinct k-mers of a sequence."""
    return {seq[i : i + K] for i in range(len(seq) - K + 1)}


def only_a(seq, start, end):
    """Whether a region holds nothing but A on this sequence."""
    return set(seq[start:end]) == {"A"}


def run_span(seq):
    """(start, end) of the poly-A run, 0-based, end exclusive."""
    start = seq.index("A")
    return start, start + seq.count("A")


def residue_axis(axis, seq):
    """One tick per residue, labelled with the residue, centred on its cell."""
    axis.set_ticks([i + 0.5 for i in range(len(seq))], list(seq), family="monospace")
    axis.set_tick_params(length=0)


def draw_sequences(ax, rows):
    """The three sequences, one residue per cell, with the poly-A run shaded."""
    for y, (label, seq) in enumerate(rows):
        start, end = run_span(seq)
        ax.add_patch(plt.Rectangle((start, y - 0.4), end - start, 0.8, color=RUN_FILL, lw=0))
        for i, residue in enumerate(seq):
            ax.text(i + 0.5, y, residue, ha="center", va="center", family="monospace")
        ax.text(len(seq) + 0.6, y, f"{len(seq)} aa", va="center")
    ax.set_xlim(0, max(len(seq) for _, seq in rows) + 3)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_yticks(range(len(rows)), [label for label, _ in rows])
    ax.tick_params(axis="y", length=0)
    ax.set_xticks([])
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    n = len(kmer_set(rows[0][1]))
    ax.set_title(
        f"All three hold the same {n} distinct {K}-mers (protein20, k = {K}); "
        "only the length of the poly-A run differs",
        loc="left",
    )


def draw_panel(ax, name):
    """Dot plot: each region is a diagonal from its first shared residue to its last."""
    target = TARGETS[name]
    t0, t1 = run_span(target)
    q0, q1 = run_span(QUERY)
    ax.axvspan(t0, t1, color=RUN_FILL, lw=0, zorder=0)
    ax.axhspan(q0, q1, color=RUN_FILL, lw=0, zorder=0)
    for qs, qe, ts, te in REGIONS[name]:
        colour = POLY_A if only_a(target, ts, te) else FLANK
        ax.plot([ts, te], [qs, qe], color=colour, solid_capstyle="butt")
    ax.set_xlim(0, len(target))
    ax.set_ylim(len(QUERY), 0)
    ax.set_aspect("equal")
    residue_axis(ax.xaxis, target)
    residue_axis(ax.yaxis, QUERY)
    ax.grid(True, color="#d9d9d9", lw=0.3)
    ax.set_axisbelow(True)
    ax.set_xlabel(f"{name} entry ({len(target)} aa)")
    ax.set_title(f"{name} entry: {len(REGIONS[name])} regions found", loc="left")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", required=True, help="path stem, without extension")
    args = parser.parse_args()

    rows = [("query", QUERY), ("short entry", TARGETS["short"]), ("long entry", TARGETS["long"])]
    assert len({frozenset(kmer_set(seq)) for _, seq in rows}) == 1, "k-mer sets differ"

    pf.use_style()
    fig = plt.figure(figsize=(pf.TWO_COLUMN_MM * pf.MM, 128 * pf.MM), layout="constrained")
    grid = fig.add_gridspec(
        4,
        2,
        height_ratios=[1.6, 0.9, 3.2, 25],
        width_ratios=[len(s) for s in TARGETS.values()],
    )
    ax_head = fig.add_subplot(grid[0, :])
    ax_key = fig.add_subplot(grid[1, :])
    ax_seq = fig.add_subplot(grid[2, :])
    axs = [fig.add_subplot(grid[3, 0])]
    axs.append(fig.add_subplot(grid[3, 1], sharey=axs[0]))
    for ax in (ax_head, ax_key):
        ax.axis("off")
    ax_head.text(0, 0.5, CONCLUSION, ha="left", va="center", fontweight="bold")

    key = [
        plt.Rectangle((0, 0), 1, 1, color=RUN_FILL),
        plt.Line2D([], [], color=FLANK),
        plt.Line2D([], [], color=POLY_A),
    ]
    ax_key.legend(
        key,
        [
            "poly-A run",
            "region that reaches past the poly-A run",
            "region of A's only (the run against itself, shifted)",
        ],
        loc="center left",
        ncol=3,
        borderaxespad=0,
    )
    draw_sequences(ax_seq, rows)
    for ax, name in zip(axs, TARGETS):
        draw_panel(ax, name)
    axs[0].set_ylabel(f"query ({len(QUERY)} aa)")
    axs[1].tick_params(axis="y", labelleft=False)
    for ax, letter in zip([ax_seq, *axs], "abc"):
        pf.panel_label(ax, letter)
    pf.save(fig, args.output)


if __name__ == "__main__":
    main()
