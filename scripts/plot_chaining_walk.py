#!/usr/bin/env python3
"""Draw why BCL2 vs CED9 at hp k=12 has one region fewer after chaining seeds per diagonal.

Positions are 0-based, as `kmerseek pair` reports them. Writes .png and .svg:

    python scripts/plot_chaining_walk.py --output docs/images/bcl2_ced9_chaining_walk
"""

import argparse

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle

DIAG_173 = "#2166ac"  # shared 12-mer on diagonal 173: BCL2 80-93 vs CED9 253-266
DIAG_89 = "#b0b0b0"  # shared 12-mer on diagonal 89: BCL2 81-93 vs CED9 170-182
HYDROPHOBIC = set("AFGILMPVWY")
BCL2 = (80, "AAGPALSPVPPVV")
CED9 = (253, "MIGAGVTAGAIGI")
MONO = dict(family="monospace", fontsize=12)


def hp(seq):
    return "".join("h" if aa in HYDROPHOBIC else "p" for aa in seq)


def draw_match(ax):
    ax.set_title("A. The match: 13 residues, two overlapping shared 12-mers", loc="left", fontsize=11)
    rows = [
        (f"BCL2 {BCL2[0]}-{BCL2[0] + 13}", BCL2[1]),
        ("BCL2 in hp", hp(BCL2[1])),
        ("CED9 in hp", hp(CED9[1])),
        (f"CED9 {CED9[0]}-{CED9[0] + 13}", CED9[1]),
    ]
    for y, (label, text) in zip([5, 4, 3, 2], rows):
        ax.text(-0.5, y, label, ha="right", va="center", fontsize=10)
        for i, ch in enumerate(text):
            ax.text(i + 0.5, y, ch, ha="center", va="center", **MONO)
    ax.text(13.3, 3.5, "13 of 13 identical in hp\n1 of 13 identical in residues", va="center", fontsize=9)
    for y, start, label in [(1.1, 0, "12-mer at BCL2 80 = CED9 253"), (0.4, 1, "12-mer at BCL2 81 = CED9 254")]:
        ax.add_patch(Rectangle((start, y - 0.2), 12, 0.4, color=DIAG_173))
        ax.text(start + 12.3, y, label, va="center", fontsize=9)
    ax.text(
        0, -0.5,
        "The 12-mer at BCL2 81 (hhhhhphhhhhh) also occurs at CED9 170-182, on another diagonal.",
        va="center", fontsize=9,
    )
    ax.set_xlim(-6, 26)
    ax.set_ylim(-1, 5.8)
    ax.axis("off")


def seed_box(ax, x, bcl2, ced9, colour):
    ax.add_patch(Rectangle((x, 0.3), 2.4, 1.0, color=colour))
    text_colour = "white" if colour == DIAG_173 else "black"
    ax.text(x + 1.2, 0.8, f"BCL2 {bcl2}\nCED9 {ced9}", ha="center", va="center", fontsize=9, color=text_colour)


def region_outline(ax, x, width):
    ax.add_patch(
        FancyBboxPatch((x - 0.15, 0.15), width + 0.3, 1.3, boxstyle="round,pad=0,rounding_size=0.2",
                       fill=False, edgecolor="black", lw=1.5)
    )


def draw_walk(ax, title, seeds, regions, notes, result):
    """seeds: (x, bcl2, ced9, colour); regions: (x, width); notes: (x, text)."""
    ax.set_title(title, loc="left", fontsize=11)
    for x, bcl2, ced9, colour in seeds:
        seed_box(ax, x, bcl2, ced9, colour)
    for x, width in regions:
        region_outline(ax, x, width)
    for x, text in notes:
        ax.annotate(text, xy=(x, 1.45), xytext=(x, 2.3), ha="center", fontsize=9,
                    arrowprops=dict(arrowstyle="-", color="black", lw=0.8))
    ax.text(15.5, 0.8, result, va="center", fontsize=10)
    ax.set_xlim(-6, 26)
    ax.set_ylim(0, 3.2)
    ax.axis("off")


def legend_handles():
    return [
        Patch(color=DIAG_173, label="shared 12-mer on diagonal 173 (BCL2 80-93 vs CED9 253-266)"),
        Patch(color=DIAG_89, label="shared 12-mer on diagonal 89 (BCL2 81-93 vs CED9 170-182)"),
        Line2D([], [], color="black", lw=1.5, label="one reported region"),
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", required=True, help="path without extension; writes .png and .svg")
    args = parser.parse_args()

    fig, axes = plt.subplots(3, 1, figsize=(10, 7.2), gridspec_kw={"height_ratios": [2.2, 1, 1], "hspace": 0.35})
    fig.suptitle("BCL2 vs CED9, hp_lehninger2, k = 12: why 14 regions became 13", fontsize=12, y=0.985)
    fig.legend(handles=legend_handles(), loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=1, fontsize=9,
               frameon=False)
    draw_match(axes[0])
    draw_walk(
        axes[1], "B. main: sorted by BCL2 position only",
        seeds=[(0, 80, 253, DIAG_173), (5.5, 81, 170, DIAG_89), (11, 81, 254, DIAG_173)],
        regions=[(0, 2.4), (5.5, 2.4), (11, 2.4)],
        notes=[(3.95, "CED9 253 to 170:\nnot +1, new region"), (9.45, "CED9 170 to 254:\nnot +1, new region")],
        result="the match is split:\ntwo regions of 12 aa",
    )
    draw_walk(
        axes[2], "C. This PR: sorted by diagonal (CED9 position minus BCL2 position), then BCL2 position",
        seeds=[(0, 81, 170, DIAG_89), (3.5, 80, 253, DIAG_173), (6.1, 81, 254, DIAG_173)],
        regions=[(0, 2.4), (3.5, 5.0)],
        notes=[(6.0, "BCL2 +1 and CED9 +1: same region")],
        result="the match is whole:\none region of 13 aa",
    )
    fig.subplots_adjust(top=0.8, bottom=0.02, left=0.02, right=0.98)
    fig.savefig(f"{args.output}.png", dpi=150)
    fig.savefig(f"{args.output}.svg")


if __name__ == "__main__":
    main()
