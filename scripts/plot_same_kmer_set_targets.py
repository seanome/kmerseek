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


def only_a(seq, start, end):
    """Whether a region holds nothing but A on this sequence."""
    return set(seq[start:end]) == {"A"}


def residue_axis(axis, seq):
    """One tick per residue, labelled with the residue, centred on its cell."""
    axis.set_ticks([i + 0.5 for i in range(len(seq))], list(seq), family="monospace")
    axis.set_tick_params(length=0)


def draw_panel(ax, name):
    """Dot plot: each region is a diagonal from its first shared residue to its last."""
    target = TARGETS[name]
    for qs, qe, ts, te in REGIONS[name]:
        colour = POLY_A if only_a(target, ts, te) else FLANK
        ax.plot([ts, te], [qs, qe], color=colour, solid_capstyle="butt")
    ax.set_xlim(0, len(target))
    ax.set_ylim(len(QUERY), 0)
    ax.set_aspect("equal")
    residue_axis(ax.xaxis, target)
    residue_axis(ax.yaxis, QUERY)
    ax.grid(True, color="#e6e6e6", lw=0.3)
    ax.set_xlabel(f"{name} entry, {len(target)} aa (poly-A run of {target.count('A')})")
    ax.set_ylabel(f"query, {len(QUERY)} aa (poly-A run of {QUERY.count('A')})")
    ax.set_title(f"{len(REGIONS[name])} regions, each cut from the {name} entry", loc="left")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", required=True, help="path stem, without extension")
    args = parser.parse_args()

    pf.use_style()
    fig, axs = pf.figure(pf.TWO_COLUMN_MM, 100, ncols=2)
    axs[0].plot([], [], color=FLANK, label="region that reaches past the poly-A run")
    axs[0].plot([], [], color=POLY_A, label="region of A's only (the run against itself, shifted)")
    for ax, name, letter in zip(axs, TARGETS, "ab"):
        draw_panel(ax, name)
        pf.panel_label(ax, letter)
    pf.shared_legend(fig, ncol=2)
    pf.save(fig, args.output)


if __name__ == "__main__":
    main()
