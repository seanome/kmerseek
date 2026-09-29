#!/usr/bin/env python3
"""Draw why the CED9-BCL2 BH1 region's E-value moves when lambda is solved from the
region's own residues instead of the two whole proteins.

    python scripts/plot_bh1_region_lambda.py -o docs/images/bh1_region_lambda.png

Reads the two fixture proteins, encodes them in the Lehninger hydrophobic/polar classes
(the `hp` alphabet), and redoes the arithmetic `test_cli_search_extend_mismatch_penalty`
pins: u, lambda, bits and E at mismatch penalty C = 2, K = 0.03, m = 280 query residues,
n = 8340 database k-mers.
"""
import argparse
import math

import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

CED9_FASTA = "tests/testdata/fasta/ced9.fasta"
BCL2_FASTA = "tests/testdata/fasta/bcl2.fasta"
# The extended BH1 region, 0-based half-open, as the search reports it.
CED9_SPAN = (162, 188)
BCL2_SPAN = (138, 164)
HYDROPHOBIC = set("AFGILMPVWY")  # Lehninger classes, as sourmash's aa_to_hp
PENALTY, K, M, N = 2.0, 0.03, 280, 8340

REGION = "#d9822b"  # the region's own residues
WHOLE = "#8c8c8c"  # the two whole proteins
H_LETTER = "#1f5fa8"  # a hydrophobic residue


def read_fasta(path):
    return "".join(line.strip() for line in open(path) if not line.startswith(">"))


def encode(seq):
    return "".join("h" if aa in HYDROPHOBIC else "p" for aa in seq)


def share_h(seq):
    return encode(seq).count("h") / len(seq)


def match_probability(q, t):
    """Chance that one position of q and one of t fall in the same class."""
    hq, ht = share_h(q), share_h(t)
    return hq * ht + (1 - hq) * (1 - ht)


def ka_lambda(u, penalty):
    """Positive root of u e^x + (1-u) e^(-penalty x) = 1; 0 past the boundary."""
    if u - penalty * (1 - u) >= 0:
        return 0.0
    f = lambda x: u * math.exp(x) + (1 - u) * math.exp(-penalty * x)
    lo, hi = 0.0, 1.0
    while f(hi) <= 1:
        hi *= 2
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if f(mid) > 1 else (mid, hi)
    return (lo + hi) / 2


def score(q, t, penalty):
    matches = sum(a == b for a, b in zip(encode(q), encode(t)))
    return matches - penalty * (len(q) - matches), matches


def draw_proteins(ax, ced9, bcl2):
    rows = [("CED9 (worm)", ced9, CED9_SPAN, 1.2), ("BCL2 (human)", bcl2, BCL2_SPAN, 0.0)]
    for name, seq, (s, e), y in rows:
        ax.plot([0, len(seq)], [y, y], color="black", lw=1.5)
        ax.add_patch(plt.Rectangle((s, y - 0.18), e - s, 0.36, color=REGION))
        ax.text(-8, y, name, ha="right", va="center", fontsize=11)
        ax.text(
            len(seq) + 6, y,
            f"{len(seq)} aa; whole protein {encode(seq).count('h')} of {len(seq)} "
            f"hydrophobic ({100 * share_h(seq):.0f}%)",
            va="center", fontsize=10, color=WHOLE,
        )
        region = seq[s:e]
        ax.text(
            (s + e) / 2, y + 0.3,
            f"region {s + 1}-{e}: {encode(region).count('h')} of {e - s} hydrophobic "
            f"({100 * share_h(region):.0f}%)",
            ha="center", va="bottom", fontsize=10, color=REGION,
        )
    ax.set_xlim(-110, 560)
    ax.set_ylim(-0.4, 1.8)
    ax.axis("off")
    ax.set_title(
        "A. The matched region (orange box on the protein line) is more hydrophobic than "
        "either protein as a whole",
        loc="left", fontsize=11,
    )


def draw_alignment(ax, q, t):
    qe, te = encode(q), encode(t)
    lines = [
        ("CED9 163", q, None),
        ("", "".join("|" if a == b else " " for a, b in zip(q, t)), None),
        ("BCL2 139", t, None),
        ("", "", None),
        ("CED9 classes", qe, None),
        ("", "".join("|" if a == b else " " for a, b in zip(qe, te)), None),
        ("BCL2 classes", te, None),
    ]
    step = 0.13
    for i, (label, text, _) in enumerate(lines):
        y = 0.95 - i * step
        ax.text(0.0, y, label, fontsize=10, va="top", family="monospace")
        for j, ch in enumerate(text):
            is_h = ch == "h" or (ch in HYDROPHOBIC and i in (0, 2))
            ax.text(0.2 + j * 0.025, y, ch, fontsize=11, va="top", family="monospace",
                    color=H_LETTER if is_h else "black")
    same_aa = sum(a == b for a, b in zip(q, t))
    same_class = sum(a == b for a, b in zip(qe, te))
    ax.text(0.88, 0.95 - 1 * step, f"{same_aa} of {len(q)} identical", fontsize=10, va="top")
    ax.text(0.88, 0.95 - 5 * step, f"{same_class} of {len(q)} same class", fontsize=10, va="top")
    s, _ = score(q, t, PENALTY)
    ax.text(0.0, 0.95 - 7.3 * step,
            f"score S = {same_class} matches - {PENALTY:g} x {len(q) - same_class} mismatches "
            f"= {s:g}. S is the same either way; only what one point of S is worth changes.",
            fontsize=10, va="top")
    ax.set_xlim(0, 1.1)
    ax.set_ylim(-0.1, 1)
    ax.axis("off")
    ax.set_title("B. The region, residues and hydrophobic (h) / polar (p) classes",
                 loc="left", fontsize=11)


def draw_lambda(ax, points):
    us = [i / 1000 for i in range(300, 667)]
    ax.plot(us, [ka_lambda(u, PENALTY) for u in us], color="black", lw=1.5)
    boundary = PENALTY / (1 + PENALTY)
    ax.axvline(boundary, color="black", ls=":", lw=1)
    ax.text(boundary - 0.005, 0.95, "u = 2/3: past here\nlambda = 0, no E-value",
            ha="right", va="top", fontsize=9)
    for label, u, lam, colour, marker in points:
        ax.plot(u, lam, marker=marker, color=colour, ms=10, ls="none")
        ax.annotate(f"{label}\nu = {u:.3f}, lambda = {lam:.3f}", (u, lam),
                    xytext=(0.36, 0.12) if colour == WHOLE else (0.565, 0.5),
                    textcoords="data", fontsize=9, color=colour,
                    arrowprops=dict(arrowstyle="-", color=colour, lw=0.8))
    ax.set_xlim(0.3, 0.72)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("u, chance two positions fall in the same class")
    ax.set_ylabel("lambda (nats per point of score)")
    ax.set_title("C. More alike by chance, so each point of score is worth less",
                 loc="left", fontsize=11)


def draw_evalue(ax, points, s):
    labels, evalues, colours = [], [], []
    for label, u, lam, colour, _ in points:
        e = K * M * N * math.exp(-lam * s)
        labels.append(f"{label}\nlambda x S = {lam:.3f} x {s:g} = {lam * s:.2f} nats")
        evalues.append(e)
        colours.append(colour)
    y = [1, 0]
    ax.barh(y, evalues, color=colours, height=0.4)
    for yi, e, label in zip(y, evalues, labels):
        ax.text(e * 1.15, yi, f"E = {e:.1f}", va="center", fontsize=10)
        ax.text(1.1, yi + 0.25, label, va="bottom", fontsize=9)
    ax.set_yticks([])
    ax.set_ylim(-0.4, 1.75)
    ax.set_xscale("log")
    ax.set_xlim(1, 1000)
    ax.set_xlabel("E-value (regions this good expected by chance), log scale")
    ax.set_title(
        f"D. E = K m n e^(-lambda S), K = {K}, m = {M}, n = {N}",
        loc="left", fontsize=11,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("-o", "--output", default="docs/images/bh1_region_lambda.png")
    args = parser.parse_args()

    ced9, bcl2 = read_fasta(CED9_FASTA), read_fasta(BCL2_FASTA)
    q, t = ced9[slice(*CED9_SPAN)], bcl2[slice(*BCL2_SPAN)]
    s, _ = score(q, t, PENALTY)
    u_whole, u_region = match_probability(ced9, bcl2), match_probability(q, t)
    points = [
        ("before: from the whole proteins", u_whole, ka_lambda(u_whole, PENALTY), WHOLE, "o"),
        ("now: from the region's own residues", u_region, ka_lambda(u_region, PENALTY),
         REGION, "s"),
    ]

    fig = plt.figure(figsize=(13, 11))
    grid = fig.add_gridspec(4, 2, height_ratios=[0.15, 0.9, 1.2, 1.6], hspace=0.45,
                            wspace=0.35)
    legend_ax = fig.add_subplot(grid[0, :])
    legend_ax.axis("off")
    legend_ax.legend(
        handles=[
            Line2D([], [], color="black", lw=1.5, label="protein (line)"),
            (Patch(color=REGION), Line2D([], [], color=REGION, marker="s", ls="none", ms=9)),
            (Patch(color=WHOLE), Line2D([], [], color=WHOLE, marker="o", ls="none", ms=9)),
            Line2D([], [], color=H_LETTER, marker="$h$", ls="none", ms=9),
        ],
        labels=[
            "protein (line)",
            "the BH1 region, and numbers from its own 26 + 26 residues (now)",
            "numbers from the two whole proteins (before)",
            "hydrophobic residue or class (blue letter)",
        ],
        handler_map={tuple: HandlerTuple(ndivide=None)},
        loc="center", ncol=2, frameon=False, fontsize=10,
    )
    fig.suptitle(
        "CED9 vs BCL2, BH1 region: same score, weaker E-value, because the region's "
        "residues match more by chance than the whole proteins do",
        fontsize=12, x=0.02, ha="left", y=0.995,
    )
    draw_proteins(fig.add_subplot(grid[1, :]), ced9, bcl2)
    draw_alignment(fig.add_subplot(grid[2, :]), q, t)
    draw_lambda(fig.add_subplot(grid[3, 0]), points)
    draw_evalue(fig.add_subplot(grid[3, 1]), points, s)
    fig.savefig(args.output, dpi=150, bbox_inches="tight")
    for label, u, lam, _, _ in points:
        print(f"{label}: u {u:.6f} lambda {lam:.6f} bits "
              f"{(lam * s - math.log(K)) / math.log(2):.4f} E {K * M * N * math.exp(-lam * s):.4f}")


if __name__ == "__main__":
    main()
