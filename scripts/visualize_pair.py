#!/usr/bin/env python3
"""Render a PNG+SVG of the k-mers one query and one target share.

Two panels, from the JSON that `kmerseek pair` writes:

1. A residue ribbon around the longest matched region: both sequences boxed residue by
   residue, each with its reduced-alphabet encoding, a tick wherever the two encodings
   agree, and the matched region outlined.
2. A dot plot of every shared k-mer, query position against target position, with every
   matched region outlined. K-mers on a region's diagonal are drawn dark, the scattered
   singles grey.

With --html, also writes a self-contained interactive page: hover a dot for the k-mer and
its position in both sequences, click a run to move the ribbon onto it, change the flank.

Usage:
    kmerseek pair --query bcl2.fasta --target ced9.fasta --ksize 12 --alphabet hp \\
        --output bcl2_vs_ced9.json
    python visualize_pair.py --pair bcl2_vs_ced9.json --output-dir pair_png/ --html
"""

import argparse
import json
import os
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from visualize_hits import INK, MUTED, SECONDARY_INK, SURFACE, safe_filename, short_label
from visualize_pair_html import render_html

# Reduced-alphabet classes get a fill colour only when there are few enough to tell apart
# at a glance (hp and hpc alphabets). Larger alphabets keep the letters and drop the fill.
CLASS_COLORS = ["#e08a1e", "#2b6cb0", "#3a9d5d", "#8e5bb5"]
MAX_COLORED_CLASSES = len(CLASS_COLORS)
RESIDUE_BOX = "#ffffff"
RESIDUE_EDGE = "#b8b7b0"

# One colour for "part of a matched region" in both panels: the dashed outline in the
# ribbon and the dot plot, and the dots that fall inside one.
REGION_COLOR = INK
SINGLE_COLOR = MUTED
# Residues of padding around a run's k-mer starts in the dot plot.
RUN_BOX_PAD = 4

# Ribbon geometry, in residue units (one residue = one x unit).
RESIDUE_INCHES = 0.2
ROW_HEIGHT = 1.0
BOX_W = 0.86
BOX_H = 0.78
FONT_PT = 8


def load_pair(path):
    with open(path) as fh:
        return json.load(fh)


def class_residues(pair):
    """{class symbol: sorted residues that map to it}, read off the two sequences rather
    than from an alphabet table, so it is right for whichever alphabet was used."""
    seen = defaultdict(set)
    for side in ("query", "target"):
        for residue, cls in zip(pair[side]["sequence"], pair[side]["encoded"]):
            seen[cls].add(residue)
    return {cls: "".join(sorted(res)) for cls, res in sorted(seen.items())}


def is_reduced(pair):
    return pair["query"]["encoded"] != pair["query"]["sequence"]


def runs(pair):
    """Matched regions made of at least two consecutive shared k-mers. `kmerseek pair` also
    reports every lone shared k-mer as a region exactly k residues long; those are drawn as
    scattered singles, not outlined."""
    return [r for r in pair["regions"] if r["length"] > pair["ksize"]]


def longest_run(pair):
    """Regions are written longest first, so the first run is the ribbon's subject."""
    return next(iter(runs(pair)), None)


def ribbon_window(region, query_len, target_len, flank):
    """Start (inclusive) and end (exclusive) of the ribbon in each sequence: the region
    plus `flank` residues either side, clipped so the window stays on the diagonal in both."""
    left = min(flank, region["query_start"], region["target_start"])
    right = min(flank, query_len - region["query_end"], target_len - region["target_end"])
    return {
        "query_start": region["query_start"] - left,
        "query_end": region["query_end"] + right,
        "target_start": region["target_start"] - left,
        "target_end": region["target_end"] + right,
    }


def window_slices(pair, window):
    q, t = pair["query"], pair["target"]
    qs, qe = window["query_start"], window["query_end"]
    ts, te = window["target_start"], window["target_end"]
    return {
        "query_seq": q["sequence"][qs:qe],
        "query_enc": q["encoded"][qs:qe],
        "target_seq": t["sequence"][ts:te],
        "target_enc": t["encoded"][ts:te],
    }


def count_agreement(a, b):
    return sum(x == y for x, y in zip(a, b))


def region_counts(pair, region):
    """Identical residues and agreeing classes inside one region."""
    q, t = pair["query"], pair["target"]
    qs, qe, ts, te = (
        region["query_start"],
        region["query_end"],
        region["target_start"],
        region["target_end"],
    )
    return {
        "identical": count_agreement(q["sequence"][qs:qe], t["sequence"][ts:te]),
        "same_class": count_agreement(q["encoded"][qs:qe], t["encoded"][ts:te]),
        "length": region["length"],
    }


def in_region(kmer, region, ksize):
    """A shared k-mer sits on a region's diagonal when its start is inside the region in
    both sequences and its diagonal offset matches."""
    q_in = region["query_start"] <= kmer["query_pos"] <= region["query_end"] - ksize
    t_in = region["target_start"] <= kmer["target_pos"] <= region["target_end"] - ksize
    same_diagonal = (kmer["target_pos"] - kmer["query_pos"]) == (
        region["target_start"] - region["query_start"]
    )
    return q_in and t_in and same_diagonal


def split_kmers_by_run(pair):
    """Shared k-mers on the diagonal of a run of consecutive shared k-mers, and the singles."""
    ksize = pair["ksize"]
    on, off = [], []
    for kmer in pair["shared_kmers"]:
        (on if any(in_region(kmer, r, ksize) for r in runs(pair)) else off).append(kmer)
    return on, off


def _text_color(hex_fill):
    r, g, b = (int(hex_fill[i : i + 2], 16) / 255 for i in (1, 3, 5))
    return "#ffffff" if 0.299 * r + 0.587 * g + 0.114 * b < 0.6 else INK


class PairPlot:
    """Draws the two panels for one pair. Kept as a class so the layout numbers and the
    colour mapping are shared between them instead of threaded through arguments."""

    def __init__(self, pair, flank=10):
        self.pair = pair
        self.flank = flank
        self.ksize = pair["ksize"]
        self.query_label = short_label(pair["query"]["name"])
        self.target_label = short_label(pair["target"]["name"])
        self.reduced = is_reduced(pair)
        self.classes = class_residues(pair) if self.reduced else {}
        self.class_fill = self._class_fill()
        self.runs = runs(pair)
        self.region = longest_run(pair)
        self.on_run, self.single = split_kmers_by_run(pair)

    def _class_fill(self):
        if len(self.classes) > MAX_COLORED_CLASSES:
            return {}
        return dict(zip(self.classes, CLASS_COLORS))

    # -- layout -------------------------------------------------------------------

    def window(self):
        if self.region is None:
            return None
        return ribbon_window(
            self.region,
            len(self.pair["query"]["sequence"]),
            len(self.pair["target"]["sequence"]),
            self.flank,
        )

    def row_labels(self, window):
        """Left-hand labels for the four ribbon rows, 1-based inclusive ranges in residues."""
        t = f"{self.target_label} {window['target_start'] + 1}–{window['target_end']}"
        q = f"{self.query_label} {window['query_start'] + 1}–{window['query_end']}"
        if not self.reduced:
            return [t, q]
        return [t, self.pair["moltype"], self.pair["moltype"], q]

    def label_margin_units(self, labels):
        """Residue units needed left of the ribbon so the longest label fits."""
        char_inches = 0.6 * FONT_PT / 72
        return max(len(s) for s in labels) * char_inches / RESIDUE_INCHES + 1.5

    # Vertical budget, top to bottom, in inches.
    TITLE_IN = 0.55
    RIBBON_LEGEND_IN = 0.3
    ROW_IN = 0.28
    REGION_LABEL_IN = 0.35
    GAP_IN = 0.1
    DOT_LEGEND_IN = 0.65
    DOT_SIDE_IN = 3.4
    XLABEL_IN = 0.65

    def n_rows(self):
        return 4 if self.reduced else 2

    def ribbon_inches(self):
        if self.region is None:
            return 0.4
        return self.RIBBON_LEGEND_IN + self.n_rows() * self.ROW_IN + self.REGION_LABEL_IN

    def figure_size(self, window):
        n = window["query_end"] - window["query_start"] if window else 40
        width = max(8.0, n * RESIDUE_INCHES + 3.0)
        height = (
            self.TITLE_IN
            + self.ribbon_inches()
            + self.GAP_IN
            + self.DOT_LEGEND_IN
            + self.DOT_SIDE_IN
            + self.XLABEL_IN
        )
        return width, height

    # -- drawing ------------------------------------------------------------------

    def draw(self):
        window = self.window()
        width, height = self.figure_size(window)
        fig = plt.figure(figsize=(width, height), facecolor=SURFACE)
        ribbon_top = height - self.TITLE_IN
        ribbon_bottom = ribbon_top - self.ribbon_inches()
        ax_ribbon = fig.add_axes(
            [0.02, ribbon_bottom / height, 0.96, self.ribbon_inches() / height]
        )
        ax_dots = fig.add_axes(
            [0.9 / width, self.XLABEL_IN / height, self.DOT_SIDE_IN / width, self.DOT_SIDE_IN / height]
        )
        for i, line in enumerate(self.title_lines()):
            fig.text(0.02, 1 - (0.12 + 0.2 * i) / height, line, ha="left", va="top", fontsize=9.5, color=INK)
        if window is None:
            self._draw_no_region(ax_ribbon)
        else:
            self._draw_ribbon(ax_ribbon, window)
        self._draw_dots(ax_dots)
        return fig

    def title_lines(self):
        first = (
            f"{self.query_label} (query) vs {self.target_label} (target), "
            f"{self.pair['moltype']} {self.ksize}-mers"
        )
        n = len(self.pair["shared_kmers"])
        if self.region is None:
            return [first, f"{n} shared, none consecutive in both sequences"]
        c = region_counts(self.pair, self.region)
        second = (
            f"{n} shared: {len(self.on_run)} in runs of consecutive k-mers, {len(self.single)} singles; "
            f"longest run {c['length']} residues, {c['identical']}/{c['length']} identical, "
            f"{c['same_class']}/{c['length']} same class"
        )
        return [first, second]

    def _draw_no_region(self, ax):
        ax.set_axis_off()
        ax.text(
            0.0,
            0.5,
            f"No run to show: no two shared {self.ksize}-mers are consecutive in both sequences.",
            transform=ax.transAxes,
            fontsize=9,
            color=SECONDARY_INK,
            va="center",
        )

    def _draw_ribbon(self, ax, window):
        s = window_slices(self.pair, window)
        n = len(s["query_seq"])
        rows = self._ribbon_rows(s)
        labels = self.row_labels(window)
        margin = self.label_margin_units(labels)
        ax.set_xlim(-margin, n + 0.5)
        # y in row units: rows sit at 0..n_rows-1, the legend above, the region label below.
        ax.set_ylim(
            -0.5 - self.REGION_LABEL_IN / self.ROW_IN,
            len(rows) - 0.5 + self.RIBBON_LEGEND_IN / self.ROW_IN,
        )
        ax.set_axis_off()
        self._draw_ribbon_legend(ax)
        for row_index, (text, encoded) in enumerate(rows):
            y = (len(rows) - 1 - row_index) * ROW_HEIGHT
            ax.text(-1.0, y, labels[row_index], ha="right", va="center", fontsize=FONT_PT, color=SECONDARY_INK)
            self._draw_row(ax, y, text, encoded)
        self._draw_agreement_ticks(ax, s, len(rows))
        self._draw_region_outline(ax, window, len(rows))

    def _ribbon_rows(self, s):
        """(text, is_encoded) per row, top to bottom: target residues, target classes,
        query classes, query residues. Without a reduced alphabet only the residue rows."""
        if not self.reduced:
            return [(s["target_seq"], False), (s["query_seq"], False)]
        return [
            (s["target_seq"], False),
            (s["target_enc"], True),
            (s["query_enc"], True),
            (s["query_seq"], False),
        ]

    def _draw_row(self, ax, y, text, encoded):
        for i, ch in enumerate(text):
            fill = self.class_fill.get(ch, RESIDUE_BOX) if encoded else RESIDUE_BOX
            edge = fill if fill != RESIDUE_BOX else RESIDUE_EDGE
            ax.add_patch(
                Rectangle((i - BOX_W / 2, y - BOX_H / 2), BOX_W, BOX_H, facecolor=fill, edgecolor=edge, linewidth=0.6)
            )
            ax.text(i, y, ch, ha="center", va="center", fontsize=FONT_PT, color=_text_color(fill) if fill != RESIDUE_BOX else INK, family="monospace")

    def _draw_agreement_ticks(self, ax, s, n_rows):
        """A tick between the two middle rows wherever the classes (or, without a reduced
        alphabet, the residues) agree."""
        top = s["target_enc"] if self.reduced else s["target_seq"]
        bottom = s["query_enc"] if self.reduced else s["query_seq"]
        upper_row = n_rows // 2  # index from the bottom of the row just above the gap
        y_top = upper_row * ROW_HEIGHT - BOX_H / 2
        y_bottom = (upper_row - 1) * ROW_HEIGHT + BOX_H / 2
        for i, (a, b) in enumerate(zip(top, bottom)):
            if a == b:
                ax.plot([i, i], [y_bottom, y_top], color=SECONDARY_INK, linewidth=0.8, solid_capstyle="butt")

    def _draw_region_outline(self, ax, window, n_rows):
        start = self.region["query_start"] - window["query_start"]
        width = self.region["length"]
        ax.add_patch(
            Rectangle(
                (start - 0.5, -BOX_H / 2 - 0.12),
                width,
                (n_rows - 1) * ROW_HEIGHT + BOX_H + 0.24,
                fill=False,
                edgecolor=REGION_COLOR,
                linestyle=(0, (4, 2)),
                linewidth=1.0,
            )
        )
        n_kmers = self.region["length"] - self.ksize + 1
        label = (
            f"{n_kmers} consecutive shared {self.ksize}-mers: {self.query_label} "
            f"{self.region['query_start'] + 1}\u2013{self.region['query_end']}, {self.target_label} "
            f"{self.region['target_start'] + 1}\u2013{self.region['target_end']}"
        )
        ax.text(start - 0.5, -BOX_H / 2 - 0.35, label, ha="left", va="top", fontsize=FONT_PT, color=REGION_COLOR)

    def _draw_ribbon_legend(self, ax):
        handles = [
            Patch(facecolor=RESIDUE_BOX, edgecolor=RESIDUE_EDGE, label="residue"),
        ]
        for cls, residues in self.classes.items():
            fill = self.class_fill.get(cls, RESIDUE_BOX)
            handles.append(Patch(facecolor=fill, edgecolor=fill if fill != RESIDUE_BOX else RESIDUE_EDGE, label=f"class {cls}: {' '.join(residues)}"))
        what = "same class" if self.reduced else "same residue"
        handles.append(Line2D([], [], color=SECONDARY_INK, linewidth=0.8, label=f"{what} in both"))
        handles.append(Line2D([], [], color=REGION_COLOR, linestyle=(0, (4, 2)), label="longest run of consecutive shared k-mers"))
        ax.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(0.0, 1.0),
            bbox_transform=ax.transAxes,
            ncol=len(handles),
            frameon=False,
            fontsize=FONT_PT,
            handlelength=1.4,
            columnspacing=1.2,
        )

    def _draw_dots(self, ax):
        ax.set_facecolor(SURFACE)
        q_len = len(self.pair["query"]["sequence"])
        t_len = len(self.pair["target"]["sequence"])
        ax.set_xlim(0, q_len)
        ax.set_ylim(0, t_len)
        ax.set_aspect("equal")
        ax.set_xlabel(f"{self.query_label} position (query, {q_len} aa)", fontsize=9)
        ax.set_ylabel(f"{self.target_label} position (target, {t_len} aa)", fontsize=9)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.tick_params(labelsize=8, colors=SECONDARY_INK)
        self._scatter(ax, self.single, SINGLE_COLOR, f"single shared {self.ksize}-mer ({len(self.single)})")
        self._scatter(ax, self.on_run, REGION_COLOR, f"shared {self.ksize}-mer in a run ({len(self.on_run)})")
        for region in self.runs:
            # Dots sit at k-mer starts, so the box surrounds those, not every residue the
            # run covers; RUN_BOX_PAD keeps a two-k-mer run's box visible around its dots.
            n_kmers = region["length"] - self.ksize + 1
            ax.add_patch(
                Rectangle(
                    (region["query_start"] - RUN_BOX_PAD, region["target_start"] - RUN_BOX_PAD),
                    n_kmers + 2 * RUN_BOX_PAD,
                    n_kmers + 2 * RUN_BOX_PAD,
                    fill=False,
                    edgecolor=REGION_COLOR,
                    linestyle=(0, (4, 2)),
                    linewidth=0.9,
                )
            )
        ax.add_patch(Rectangle((0, 0), 0, 0, fill=False, edgecolor=REGION_COLOR, linestyle=(0, (4, 2)), label=f"run of consecutive shared k-mers ({len(self.runs)})"))
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), frameon=False, fontsize=8, ncol=1, borderaxespad=0)

    def _scatter(self, ax, kmers, color, label):
        # A k-mer covers k residues; the dot sits on its start, plus half a residue so a
        # k-mer starting at 0-based position 0 is drawn inside the axes.
        xs = [k["query_pos"] + 0.5 for k in kmers]
        ys = [k["target_pos"] + 0.5 for k in kmers]
        ax.scatter(xs, ys, s=9, color=color, label=label, linewidths=0, zorder=3)


def plot_pair(pair, output_paths, flank=10, dpi=200):
    fig = PairPlot(pair, flank=flank).draw()
    for path in output_paths:
        fig.savefig(path, dpi=dpi, facecolor=SURFACE)
    plt.close(fig)


def write_html(pair, path, flank=10):
    title = f"{short_label(pair['query']['name'])} vs {short_label(pair['target']['name'])} shared k-mers"
    with open(path, "w") as fh:
        fh.write(render_html(pair, title, flank=flank))


def output_basename(pair):
    return f"{safe_filename(pair['query']['name'])}_vs_{safe_filename(pair['target']['name'])}.{pair['moltype']}.k{pair['ksize']}"


def _build_arg_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pair", required=True, help="JSON written by `kmerseek pair`")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--flank", type=int, default=10, help="residues shown either side of the longest matched region (default 10)")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument(
        "--html",
        action="store_true",
        help="also write a self-contained interactive HTML page (hover a dot for the k-mer, click a run to show it)",
    )
    return p


def main():
    args = _build_arg_parser().parse_args()
    pair = load_pair(args.pair)
    os.makedirs(args.output_dir, exist_ok=True)
    base = os.path.join(args.output_dir, output_basename(pair))
    plot_pair(pair, [base + ".png", base + ".svg"], flank=args.flank, dpi=args.dpi)
    print(base + ".png")
    print(base + ".svg")
    if args.html:
        write_html(pair, base + ".html", flank=args.flank)
        print(base + ".html")


if __name__ == "__main__":
    main()
