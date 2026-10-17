#!/usr/bin/env python3
"""Draw topi pancreatic ribonuclease (P00659) against itself at protein20 k=10, before and
after `find_matched_regions` keeps each shared (query, target) position once.

A 10-mer holding n residues that are B or Z is sketched under all 2^n readings, and both
copies hold every reading, so its start position is listed 2^n times on the main diagonal.
The run walk joins two neighbours only when both positions go up by one, so before the
fix it broke at every repeat. The walk below is the one in `find_matched_regions`; its
output matches `cargo test ... an_ambiguous_residue_shared_with_itself` on both builds
(427 regions before, 1 after). Positions are 0-based, as `kmerseek pair` reports them.

    python scripts/plot_ambiguous_position_dedup.py \
        --output docs/images/topi_self_ambiguous_position_dedup
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import pubfig as pf  # noqa: E402

TOPI = (
    "KESAAAKFZRZHMBSSTSSASSSBYCBZMMKSRNLTQDRCKPVBTFVHZSLABVZAVCSZ"
    "KBVACKBGZTBCYZSYSTMSITBCRZTGSSKYPBCAYKTTQAKKHIIVACZGBPYVPVHF"
    "BASV"
)
K = 10
AMBIGUOUS = "BZJ"
REGION = pf.OKABE_ITO["blue"]
AMBIGUOUS_COLOR = pf.OKABE_ITO["vermillion"]


def listed_starts(seq, k):
    """Each 10-mer start once per reading: 2^n times for n ambiguous residues in it."""
    starts = []
    for q in range(len(seq) - k + 1):
        starts += [q] * 2 ** sum(aa in AMBIGUOUS for aa in seq[q : q + k])
    return starts


def walk(starts, k):
    """Regions as (start, end): runs whose neighbours both go up by one."""
    regions, i = [], 0
    while i < len(starts):
        j = i + 1
        while j < len(starts) and starts[j] == starts[j - 1] + 1:
            j += 1
        regions.append((starts[i], starts[i] + (j - i) + k - 1))
        i = j
    return regions


def pileup_rows(regions):
    """Put each region on the lowest row where it overlaps nothing, as in a read pileup."""
    row_ends, rows = [], []
    for start, end in regions:
        row = next((r for r, e in enumerate(row_ends) if e <= start), len(row_ends))
        if row == len(row_ends):
            row_ends.append(end)
        row_ends[row] = end
        rows.append(row)
    return rows


def draw_protein(ax):
    ax.plot([0, len(TOPI)], [1, 1], color="black", lw=0.75)
    ambiguous = [i for i, aa in enumerate(TOPI) if aa in AMBIGUOUS]
    ax.plot([i + 0.5 for i in ambiguous], [1.6] * len(ambiguous), "v",
            color=AMBIGUOUS_COLOR, ms=3, label="B or Z (ambiguous residue)")
    for i, aa in enumerate(TOPI):
        ax.text(i + 0.5, 0.2, aa, family="monospace", fontsize=5, ha="center", va="center")
    ax.set_ylim(-0.4, 2.1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_title(f"Topi ribonuclease (P00659), {len(TOPI)} aa: "
                 f"{TOPI.count('B')} B and {TOPI.count('Z')} Z", loc="left")


def draw_listed(ax, starts):
    counts = {q: starts.count(q) for q in sorted(set(starts))}
    ax.plot([q + 0.5 for q in counts], list(counts.values()), "o", color="black", ms=2,
            label="Times a 10-mer start is listed")
    ax.set_yticks([1, 2, 4, 8, 16])
    ax.set_ylabel("Times listed")
    ax.set_title(f"Each 10-mer start is listed $2^n$ times, n = number of B or Z in it: "
                 f"{len(starts)} entries for {len(counts)} starts", loc="left")


def draw_regions(ax, regions, title):
    rows = pileup_rows(regions)
    for (start, end), row in zip(regions, rows):
        ax.plot([start, end], [row, row], color=REGION, lw=1.2, solid_capstyle="butt",
                label="Matched region reported")
    ax.set_ylim(max(rows) + 1, -1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_title(title, loc="left")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", required=True, help="path stem; writes .svg .png")
    stem = parser.parse_args().output

    starts = listed_starts(TOPI, K)
    before, after = walk(starts, K), walk(sorted(set(starts)), K)
    single = sum(end - start == K for start, end in before)
    assert (len(before), single, after) == (427, 326, [(0, len(TOPI))]), (len(before), single)

    pf.use_style()
    fig, axs = pf.figure(pf.TWO_COLUMN_MM, 150, nrows=4, sharex=True,
                         height_ratios=[1.2, 2, 6, 0.6])
    draw_protein(axs[0])
    draw_listed(axs[1], starts)
    draw_regions(axs[2], before, f"Before this PR: {len(before)} regions, "
                 f"{single} of them a single 10-mer long (drawn so no two overlap on a row)")
    draw_regions(axs[3], after, f"After this PR: {len(after)} region, residues 0-{len(TOPI)}, "
                 "the right answer for a protein against itself")
    axs[3].set_xlabel("Topi residue, query and target (aa, 0-based); kmerseek protein20, k=10")
    axs[3].set_xlim(0, len(TOPI))
    pf.shared_legend(fig)
    for ax, letter in zip(axs, "abcd"):
        pf.panel_label(ax, letter)
    pf.save(fig, stem, formats=("svg", "png"))


if __name__ == "__main__":
    main()
