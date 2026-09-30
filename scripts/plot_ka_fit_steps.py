#!/usr/bin/env python3
"""Draw how `kmerseek search` fits r_database and K before a search: the steps as boxes
(top panel), and one real fit (bottom panel) from a survival CSV written by
`kmerseek index --ka-survival-out` with `--ka-null database`.

    python scripts/plot_ka_fit_steps.py docs/images/ka_fit_curves/scope40_database.csv \\
        -o docs/images/ka_fit_steps.png --subtitle "SCOPe40, hp_thomas_dill2 k = 12, C = 2, X = 8"
"""

import argparse
import csv
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

REAL = "#1f77b4"
REFERENCE = "#7f7f7f"
FITTED = "#d62728"
WINDOW = "#e0e0e0"

REAL_ROW = [
    "200 database\nsequences\n(--ka-null database)",
    "search each one\nagainst the index;\ndrop its hit on itself",
    "score every region:\nx = λ_pair × S\n(nats)",
    "count regions in\neach 0.5-nat bin of x",
]
REFERENCE_ROW = [
    "the same 200,\nshuffled keeping\ndipeptides\n(--ka-reference)",
    "search each one\nagainst the index",
    "score every region:\nx = λ_pair × S",
    "count regions in\neach bin of x",
]
FIT_BOX = "fit a line to ln(count)\nagainst x, below the bin\nwhere real counts rise\nabove shuffled ones.\n−slope = r_database\nheight → K"


def box(ax, x, y, text, color):
    ax.add_patch(FancyBboxPatch((x, y), 1.7, 0.8, boxstyle="round,pad=0.04", fc="white", ec=color, lw=1.5))
    ax.text(x + 0.85, y + 0.4, text, ha="center", va="center", fontsize=9)


def arrow(ax, start, end, color):
    ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", color=color, lw=1.3))


def draw_steps(ax):
    for row, (texts, color, y) in enumerate([(REAL_ROW, REAL, 1.2), (REFERENCE_ROW, REFERENCE, 0.0)]):
        for i, text in enumerate(texts):
            box(ax, i * 2.1, y, text, color)
            if i:
                arrow(ax, (i * 2.1 - 0.4, y + 0.4), (i * 2.1, y + 0.4), color)
        arrow(ax, (3 * 2.1 + 1.7, y + 0.4), (8.7, 1.0), color)
    ax.add_patch(FancyBboxPatch((8.7, 0.35), 2.2, 1.3, boxstyle="round,pad=0.04", fc="white", ec=FITTED, lw=1.5))
    ax.text(9.8, 1.0, FIT_BOX, ha="center", va="center", fontsize=9)
    ax.set_xlim(-0.1, 11.0)
    ax.set_ylim(-0.1, 2.1)
    ax.axis("off")


def load(path):
    rows = list(csv.DictReader(open(path)))
    x = [float(r["x"]) for r in rows]

    def per_bin(col):
        at_least = [int(r[col] or 0) for r in rows]
        return [a - (at_least[i + 1] if i + 1 < len(at_least) else 0) for i, a in enumerate(at_least)]

    meta = rows[0]
    window = [float(r["x"]) for r in rows if r["in_fit"] == "true"]
    return x, per_bin("n_regions_at_least"), per_bin("reference_n_regions_at_least"), window, meta


def draw_fit(ax, path):
    x, real, reference, window, meta = load(path)
    r, k, width = float(meta["r_database"]), float(meta["k"]), float(meta["bin_width"])
    scale = k * float(meta["query_residues"]) * float(meta["database_kmers"]) * (1 - math.exp(-r * width))
    ax.axvspan(window[0], window[-1] + width, color=WINDOW, zorder=0)
    ax.plot([s for s, c in zip(x, reference) if c > 0], [c for c in reference if c > 0], "^", mfc="none", color=REFERENCE, ms=6, zorder=2)
    ax.plot(x, real, "o", color=REAL, ms=4, zorder=3)
    line_x = [s for s in x if s >= 5]
    ax.plot(line_x, [scale * math.exp(-r * s) for s in line_x], color=FITTED, lw=1.5, zorder=4)
    ax.set_yscale("log")
    ax.set_ylim(0.5, 1e8)
    ax.set_xlabel("x = λ_pair × S, the region's score in nats")
    ax.set_ylabel("regions in the 0.5-nat bin")
    stop = window[-1] + width
    ax.annotate(
        "the fit stops here: above it the\nreal counts rise over the shuffled\nones because of relatives",
        xy=(stop, 3e3), xytext=(stop + 6, 1e5), fontsize=9,
        arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
    )
    ax.spines[["top", "right"]].set_visible(False)
    return r, k


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--subtitle", default="")
    args = parser.parse_args()
    fig = plt.figure(figsize=(11, 8.2))
    legend_ax = fig.add_axes([0.06, 0.9, 0.9, 0.06])
    steps_ax = fig.add_axes([0.02, 0.56, 0.96, 0.32])
    fit_ax = fig.add_axes([0.09, 0.08, 0.86, 0.4])
    draw_steps(steps_ax)
    r, k = draw_fit(fit_ax, args.csv)
    handles = [
        plt.Line2D([], [], marker="o", ls="none", color=REAL, label="database sequences, as they are (blue)"),
        plt.Line2D([], [], marker="^", ls="none", mfc="none", color=REFERENCE, label="the same sequences shuffled keeping dipeptides (grey)"),
        plt.Line2D([], [], color=FITTED, label=f"fitted line: r_database = {r:.3f}, K = {k:.4f}"),
        plt.Rectangle((0, 0), 1, 1, color=WINDOW, label="bins the line is fitted to"),
    ]
    legend_ax.legend(handles=handles, loc="center", ncol=2, frameon=False, fontsize=9)
    legend_ax.axis("off")
    fig.suptitle(f"How r_database and K are fitted before a search. {args.subtitle}", fontsize=11, y=0.99)
    fig.savefig(args.output, dpi=150)


if __name__ == "__main__":
    main()
