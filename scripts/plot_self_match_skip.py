#!/usr/bin/env python3
"""Draw which hits kmerseek reports for an index holding two proteins, BCL-2 (BCL2_HUMAN)
and Ced-9 (CED9_CAEEL), at hp_lehninger2 k=12.

Panel a: BCL-2 from a FASTA file searched against the index. The query has the same
sequence as the BCL2_HUMAN entry, so their keys are equal. Before this change, compare()
dropped any target with the query's key, so BCL2_HUMAN was missing. Now the entry is
reported.

Panel b: the index searched against itself (all-vs-all). Each entry is also a query, and its
hit on itself is still dropped, leaving two rows.

The hits are the ones
`cargo test --no-default-features --lib query_identical_to_an_entry_hits_it_unless_all_vs_all`
asserts.

    python scripts/plot_self_match_skip.py \
        --output docs/images/self_match_skip_bcl2_ced9_hp_lehninger2_k12
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import pubfig as pf  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

TARGETS = ["BCL2_HUMAN", "CED9_CAEEL"]
REPORTED = pf.OKABE_ITO["blue"]
COLUMN_STEP = 2  # wide enough that the entry names under the squares do not overlap

# (row label, reported per target) for each panel.
FASTA_QUERY = [
    ("before: BCL-2 query", [False, True]),
    ("after: BCL-2 query", [True, True]),
]
ALL_VS_ALL = [
    ("BCL2_HUMAN as query", [False, True]),
    ("CED9_CAEEL as query", [True, False]),
]


def draw_grid(ax, rows, title):
    """One square per (row, target): filled when reported, open when dropped."""
    for y, (_, reported) in enumerate(rows):
        for x, hit in enumerate(reported):
            ax.add_patch(
                Rectangle(
                    (x * COLUMN_STEP - 0.4, y - 0.4),
                    0.8,
                    0.8,
                    facecolor=REPORTED if hit else "none",
                    edgecolor=REPORTED if hit else pf.GREY,
                )
            )
    ax.set_xlim(-1, (len(TARGETS) - 1) * COLUMN_STEP + 1)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xticks([i * COLUMN_STEP for i in range(len(TARGETS))], TARGETS)
    ax.set_yticks(range(len(rows)), [label for label, _ in rows])
    ax.set_xlabel("database entry")
    ax.set_title(title)
    ax.set_aspect("equal")
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", required=True, help="path stem, no extension")
    args = parser.parse_args()

    pf.use_style()
    fig, (a, b) = pf.figure(pf.TWO_COLUMN_MM, 45, ncols=2)
    draw_grid(a, FASTA_QUERY, "FASTA query against the index")
    draw_grid(b, ALL_VS_ALL, "index against itself (all-vs-all)")
    a.plot([], [], "s", color=REPORTED, label="hit reported")
    a.plot([], [], "s", markerfacecolor="none", color=pf.GREY, label="hit dropped")
    pf.shared_legend(fig)
    for ax, letter in zip((a, b), "ab"):
        pf.panel_label(ax, letter)
    pf.save(fig, args.output)


if __name__ == "__main__":
    main()
