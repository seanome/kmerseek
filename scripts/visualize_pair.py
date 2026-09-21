#!/usr/bin/env python3
"""Render a PNG+SVG of the k-mers one query and one target share.

From the JSON that `kmerseek pair` writes:

1. A dot plot of every shared k-mer, query position against target position. Each
   protein is drawn as a line with boxes for its domains along its axis (with --domains),
   and each domain's span is shaded across the plot, so a run sits in a named cell such
   as "Bcl-2 x Bcl-2" without reading coordinates. Runs of two or more consecutive shared
   k-mers are numbered diagonal segments; lone shared k-mers are dots.
2. One alignment block per run, longest first, numbered to match the dot plot: the query
   row, the identical residues written between the rows, and the target row, with
   1-based coordinates at both ends and residue boxes coloured by alphabet class. The
   header gives the two regions, the length, the identical residues and, for a
   hydrophobic/polar alphabet, how many residues are polar (a low-complexity flag).

--domains takes Pfam-style tables (TSV, CSV or parquet) with one row per domain: a
protein column (accession, protein or name), start and end columns (domain_start/
domain_end or start/end, 1-based inclusive), and a name column (name, pfam_name or
pfam_id). Proteins match by full FASTA header, first token, or UniProt accession.
--html also writes a self-contained interactive page.

Usage:
    kmerseek pair --query bcl2.fasta --target ced9.fasta --ksize 12 --alphabet hp \\
        --output bcl2_vs_ced9.json
    python visualize_pair.py --pair bcl2_vs_ced9.json --output-dir pair_png/ \\
        --domains pfam_domains.tsv --html
"""

import argparse
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pair_model import build_model, load_domains, run_header, title_lines
from structure_alignment import align_pair, find_aligner
from visualize_hits import INK, SECONDARY_INK, SURFACE, safe_filename
from hits_page import render_page

# Residue box fill and edge per alphabet class, in class-symbol order. Only alphabets with
# this few classes get colour; larger alphabets keep the letters and plain boxes.
CLASS_STYLES = [("#f3e3c3", "#c9a15a"), ("#dde9f8", "#7fa6d6"), ("#dff0e2", "#6fae7c"), ("#ece0f3", "#a98bc4")]
PLAIN_BOX = ("#ffffff", "#b8b7b0")

DOMAIN_FILL, DOMAIN_EDGE = "#e8e7e2", "#8d8c86"
SPAN_SHADE = "#e8e7e2"
RUN_COLOR = INK
SINGLE_COLOR = "#6e6d68"
# USalign's residue pairs, drawn as a thin line under the runs.
STRUCTURE_COLOR = "#c0392b"

WRAP = 60  # residues per alignment line before wrapping
RESIDUE_IN = 0.155
ROW_IN = 0.2
PLOT_IN = 3.6
TRACK_IN = 0.42
FONT = 8


def path_segments(pairs):
    """USalign's residue pairs as runs of consecutive pairs, each a list of query positions
    and the matching target positions, so gaps break the drawn line."""
    segments, current = [], []
    for q, t, _ in pairs:
        if current and (q, t) != (current[-1][0] + 1, current[-1][1] + 1):
            segments.append(current)
            current = []
        current.append((q, t))
    if current:
        segments.append(current)
    return [([q for q, _ in seg], [t for _, t in seg]) for seg in segments]


def block_rows(block):
    """The rows an alignment block shows: the exact run with any flank, and which of its
    columns are the run."""
    left = block["query_start"] - block["window"]["query_start"]
    return {k: block[k] for k in ("query_row", "target_row", "query_enc", "target_enc", "middle")} | {
        "query_start": block["window"]["query_start"],
        "target_start": block["window"]["target_start"],
        "run_columns": [left, left + block["length"]],
    }


def axis_ticks(length):
    """1, every 100 residues (200 past 600 aa, 500 past 1500), and the length; a round tick
    within 5% of the length would overlap its label."""
    step = 100 if length <= 600 else 200 if length <= 1500 else 500
    rounds = [t for t in range(step, length, step) if length - t > length * 0.05]
    return [1, *rounds, length]


def load_pair(path):
    with open(path) as fh:
        return json.load(fh)


class PairFigure:
    """Draws the dot plot and the stacked alignment blocks from one figure model."""

    def __init__(self, model):
        self.m = model
        self.styles = self._styles()

    def _styles(self):
        classes = [c["symbol"] for c in self.m["classes"]]
        if 0 < len(classes) <= len(CLASS_STYLES):
            return dict(zip(classes, CLASS_STYLES))
        return {}

    def style(self, cls):
        return self.styles.get(cls, PLAIN_BOX)

    # -- layout in inches, top to bottom --

    def has_domains(self):
        return bool(self.m["query"]["domains"] or self.m["target"]["domains"])

    def legend_rows(self):
        return len(self._legend_handles())

    @staticmethod
    def chunks(block):
        return -(-len(block_rows(block)["query_row"]) // WRAP)

    def height(self):
        h = 0.15 + 0.22 * len(title_lines(self.m)) + 0.2 * self.legend_rows() + 0.35
        h += TRACK_IN + 0.15 + PLOT_IN + 0.55
        for b in self.m["runs"]:
            h += 0.32 + self.chunks(b) * (3 * ROW_IN + 0.2)
        return h + 0.2

    def width(self):
        """Wide enough for the title (about 0.6 * 9.5 pt per character), the widest
        alignment line, and the dot plot with its target track."""
        longest = max((len(block_rows(b)["query_row"]) for b in self.m["runs"]), default=0)
        title_in = max(len(line) for line in title_lines(self.m)) * 0.6 * 9.5 / 72 + 0.3
        return max(7.5, title_in, 1.3 + min(longest, WRAP) * RESIDUE_IN + 1.0)

    def draw(self):
        w, h = self.width(), self.height()
        self.w, self.h = w, h
        fig = plt.figure(figsize=(w, h), facecolor=SURFACE)
        y = h - 0.1
        y = self._draw_title(fig, y)
        y = self._draw_legend(fig, y)
        y = self._draw_dotplot(fig, y)
        for block in self.m["runs"]:
            y = self._draw_block(fig, y, block)
        return fig

    def _text(self, fig, x_in, y_in, text, **kw):
        fig.text(x_in / self.w, y_in / self.h, text, **kw)

    def _axes(self, fig, x_in, y_in, w_in, h_in):
        return fig.add_axes([x_in / self.w, y_in / self.h, w_in / self.w, h_in / self.h])

    # -- title and legend --

    def _draw_title(self, fig, y):
        lines = title_lines(self.m)
        self._text(fig, 0.1, y, lines[0], fontsize=9.5, color=INK, va="top")
        for i, line in enumerate(lines[1:], start=1):
            self._text(fig, 0.1, y - 0.22 * i, line, fontsize=8.5, color=SECONDARY_INK, va="top")
        return y - 0.22 * len(lines) - 0.12

    def _legend_handles(self):
        k = self.m["ksize"]
        handles = [Patch(facecolor=self.style(c["symbol"])[0], edgecolor=self.style(c["symbol"])[1], label=c["label"]) for c in self.m["classes"]]
        handles = handles or [Patch(facecolor=PLAIN_BOX[0], edgecolor=PLAIN_BOX[1], label="residue")]
        handles.append(Line2D([], [], color=RUN_COLOR, linewidth=2.2, label=f"run of 2 or more consecutive shared {k}-mers ({len(self.m['runs'])}), numbered; underlined in its alignment"))
        handles.append(Line2D([], [], color=SINGLE_COLOR, marker="o", linestyle="none", markersize=4, label=f"single shared {k}-mer ({len(self.m['singles'])})"))
        if self.has_domains():
            handles.append(Patch(facecolor=DOMAIN_FILL, edgecolor=DOMAIN_EDGE, label="protein, with its domains as boxes; each domain's span shaded across the plot"))
        if self.m.get("structure"):
            st = self.m["structure"]
            handles.append(Line2D([], [], color=STRUCTURE_COLOR, linewidth=1.2, label=f"residue pairs from superposing the two AlphaFold models with {st['aligner']} (TM-score {st['tm_score_query']:.2f})"))
        handles.append(Line2D([], [], color=INK, marker="$\\mathtt{G}$", linestyle="none", markersize=6, label="identical residue, written between the rows"))
        return handles

    def _draw_legend(self, fig, y):
        rows = self.legend_rows()
        ax = self._axes(fig, 0.1, y - 0.2 * rows, self.w - 0.2, 0.2 * rows)
        ax.set_axis_off()
        ax.legend(handles=self._legend_handles(), loc="upper left", ncol=1, frameon=False, fontsize=FONT, handlelength=1.6, borderaxespad=0, labelspacing=0.35)
        return y - 0.2 * rows - 0.25

    # -- dot plot with protein tracks --

    def _draw_dotplot(self, fig, y):
        left = 0.85
        q, t = self.m["query"], self.m["target"]
        self._draw_track(self._axes(fig, left, y - TRACK_IN, PLOT_IN, TRACK_IN), q, horizontal=True)
        top = y - TRACK_IN - 0.05
        ax = self._axes(fig, left, top - PLOT_IN, PLOT_IN, PLOT_IN)
        self._draw_track(self._axes(fig, left + PLOT_IN + 0.05, top - PLOT_IN, TRACK_IN, PLOT_IN), t, horizontal=False)
        self._draw_spans(ax)
        self._draw_marks(ax)
        self._style_axes(ax, q, t)
        return top - PLOT_IN - 0.55

    def _style_axes(self, ax, q, t):
        ax.set_xlim(1, q["length"])
        ax.set_ylim(1, t["length"])
        ax.set_xticks(axis_ticks(q["length"]))
        ax.set_yticks(axis_ticks(t["length"]))
        ax.tick_params(labelsize=FONT, colors=SECONDARY_INK)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.set_xlabel(f"{q['label']} position (aa)", fontsize=FONT + 1)
        ax.set_ylabel(f"{t['label']} position (aa)", fontsize=FONT + 1)

    def _draw_track(self, ax, side, horizontal):
        """The protein as a line with a box per domain, labelled."""
        ax.set_axis_off()
        L = side["length"]
        if horizontal:
            ax.set_xlim(1, L)
            ax.set_ylim(0, 1)
            ax.plot([1, L], [0.3, 0.3], color=DOMAIN_EDGE, linewidth=1)
        else:
            ax.set_ylim(1, L)
            ax.set_xlim(0, 1)
            ax.plot([0.3, 0.3], [1, L], color=DOMAIN_EDGE, linewidth=1)
        stagger = self._labels_collide(side) if horizontal else False
        for i, d in enumerate(side["domains"]):
            self._draw_domain(ax, d, horizontal, lift=0.3 * (i % 2) if stagger else 0)

    def _draw_domain(self, ax, d, horizontal, lift):
        span = d["end"] - d["start"] + 1
        mid = (d["start"] + d["end"]) / 2
        if horizontal:
            ax.add_patch(Rectangle((d["start"], 0.1), span, 0.4, facecolor=DOMAIN_FILL, edgecolor=DOMAIN_EDGE, linewidth=0.8))
            ax.text(mid, 0.6 + lift, d["name"], ha="center", va="bottom", fontsize=FONT, color=SECONDARY_INK)
        else:
            ax.add_patch(Rectangle((0.1, d["start"]), 0.4, span, facecolor=DOMAIN_FILL, edgecolor=DOMAIN_EDGE, linewidth=0.8))
            ax.text(0.6, mid, d["name"], ha="left", va="center", fontsize=FONT, color=SECONDARY_INK)

    def _labels_collide(self, side):
        """Whether two neighbouring domain labels on the horizontal track would overlap,
        at about 0.6 * FONT points per character."""
        aa_per_pt = side["length"] / (PLOT_IN * 72)
        centres = [((d["start"] + d["end"]) / 2, len(d["name"])) for d in side["domains"]]
        for (c1, n1), (c2, n2) in zip(centres, centres[1:]):
            if (c2 - c1) < ((n1 + n2) / 2 * 0.6 * FONT + 4) * aa_per_pt:
                return True
        return False

    def _draw_spans(self, ax):
        for d in self.m["query"]["domains"]:
            ax.axvspan(d["start"], d["end"] + 1, facecolor=SPAN_SHADE, alpha=0.6, linewidth=0, zorder=0)
        for d in self.m["target"]["domains"]:
            ax.axhspan(d["start"], d["end"] + 1, facecolor=SPAN_SHADE, alpha=0.6, linewidth=0, zorder=0)

    def _draw_structure_path(self, ax):
        if not self.m.get("structure"):
            return
        for qs, ts in path_segments(self.m["structure"]["pairs"]):
            ax.plot([p + 1 for p in qs], [p + 1 for p in ts], color=STRUCTURE_COLOR, linewidth=1.2, solid_capstyle="round", zorder=2)

    def _draw_marks(self, ax):
        self._draw_structure_path(ax)
        k = self.m["ksize"]
        # A single k-mer is a dot at the centre of the k residues it covers.
        xs = [s["query_pos"] + (k + 1) / 2 for s in self.m["singles"]]
        ys = [s["target_pos"] + (k + 1) / 2 for s in self.m["singles"]]
        ax.scatter(xs, ys, s=12, color=SINGLE_COLOR, linewidths=0, zorder=3)
        for b in self.m["runs"]:
            x = (b["query_start"] + 1, b["query_end"])
            yy = (b["target_start"] + 1, b["target_end"])
            ax.plot(x, yy, color=RUN_COLOR, linewidth=2.2, solid_capstyle="round", zorder=4)
            ax.annotate(str(b["number"]), (x[1], yy[1]), xytext=(4, 2), textcoords="offset points", fontsize=FONT + 2, fontweight="bold", color=RUN_COLOR)

    # -- alignment blocks --

    def _draw_block(self, fig, y, block):
        self._text(fig, 0.1, y - 0.05, run_header(block), fontsize=FONT + 2, color=INK, va="top")
        y -= 0.32
        rows = block_rows(block)
        for start in range(0, len(rows["query_row"]), WRAP):
            y = self._draw_chunk(fig, y, rows, start)
        return y

    def _draw_chunk(self, fig, y, rows, start):
        sl = slice(start, start + WRAP)
        n = len(rows["query_row"][sl])
        ax = self._axes(fig, 1.3, y - 3.4 * ROW_IN, n * RESIDUE_IN, 3.4 * ROW_IN)
        ax.set_axis_off()
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_ylim(-0.9, 2.5)
        self._draw_run_bar(ax, rows["run_columns"], start, n)
        self._draw_row(ax, 2, rows["query_row"][sl], rows["query_enc"][sl])
        self._draw_row(ax, 0, rows["target_row"][sl], rows["target_enc"][sl])
        for i, ch in enumerate(rows["middle"][sl]):
            if ch != " ":
                ax.text(i, 1, ch, ha="center", va="center", fontsize=FONT, color=INK, family="monospace")
        self._draw_coordinates(ax, rows, start, n)
        return y - 3.4 * ROW_IN - 0.12

    def _draw_run_bar(self, ax, run_columns, start, n):
        """The run's columns, marked under the target row with the dot plot's run bar."""
        first, last = max(run_columns[0], start) - start, min(run_columns[1], start + n) - start
        if last > first:
            ax.plot([first - 0.4, last - 0.6], [-0.7, -0.7], color=RUN_COLOR, linewidth=2.2, solid_capstyle="butt")

    def _draw_coordinates(self, ax, rows, start, n):
        """1-based first and last residue of each row in this chunk."""
        for y, label, origin in (
            (2, self.m["query"]["label"], rows["query_start"]),
            (0, self.m["target"]["label"], rows["target_start"]),
        ):
            first, last = origin + start, origin + start + n
            ax.text(-0.9, y, f"{label} {first + 1}", ha="right", va="center", fontsize=FONT, color=SECONDARY_INK)
            ax.text(n - 0.1, y, str(last), ha="left", va="center", fontsize=FONT, color=SECONDARY_INK)

    def _draw_row(self, ax, y, text, encoded):
        for i, (ch, cls) in enumerate(zip(text, encoded)):
            fill, edge = self.style(cls)
            ax.add_patch(Rectangle((i - 0.42, y - 0.4), 0.84, 0.8, facecolor=fill, edgecolor=edge, linewidth=0.6))
            ax.text(i, y, ch, ha="center", va="center", fontsize=FONT, color=INK, family="monospace")


def plot_pair(model, output_paths, dpi=200):
    fig = PairFigure(model).draw()
    for path in output_paths:
        fig.savefig(path, dpi=dpi, facecolor=SURFACE)
    plt.close(fig)


def write_html(model, path):
    """The pair as a one-row page of the search report template, the row open."""
    title = f"{model['query']['label']} vs {model['target']['label']} kmerseek pair"
    with open(path, "w") as fh:
        fh.write(render_page(title, {"pair": model}))


def structure_for(args, pair):
    """The USalign superposition for the pair when --structures is given, else None. Missing
    files or aligner are reported, not fatal."""
    if not args.structures:
        return None
    aligner = find_aligner(args.aligner)
    if aligner is None:
        print("no USalign or TMalign found; skipping the superposition", file=sys.stderr)
        return None
    structure = align_pair(aligner, args.structures, pair["query"]["name"], pair["target"]["name"])
    if structure is None:
        print(f"no structure file for both proteins in {args.structures}; skipping the superposition", file=sys.stderr)
    return structure


def output_basename(pair):
    return f"{safe_filename(pair['query']['name'])}_vs_{safe_filename(pair['target']['name'])}.{pair['moltype']}.k{pair['ksize']}"


def _build_arg_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pair", required=True, help="JSON written by `kmerseek pair`")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--domains", nargs="*", default=[], metavar="TABLE", help="domain tables (TSV, CSV or parquet) for either protein; see the module docstring for columns")
    p.add_argument("--flank", type=int, default=0, help="residues shown either side of each run (default 0); with a flank the middle line uses `:` for same class")
    p.add_argument("--structures", metavar="DIR", help="directory of AlphaFold or PDB files; with USalign or TM-align, their residue pairs are drawn across the dot plot")
    p.add_argument("--aligner", help="USalign or TMalign binary (default: found on PATH)")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--html", action="store_true", help="also write a self-contained interactive HTML page")
    return p


def main():
    args = _build_arg_parser().parse_args()
    pair = load_pair(args.pair)
    model = build_model(
        pair,
        load_domains(args.domains),
        flank=args.flank,
        structure=structure_for(args, pair),
    )
    os.makedirs(args.output_dir, exist_ok=True)
    base = os.path.join(args.output_dir, output_basename(pair))
    plot_pair(model, [base + ".png", base + ".svg"], dpi=args.dpi)
    print(base + ".png")
    print(base + ".svg")
    if args.html:
        write_html(model, base + ".html")
        print(base + ".html")


if __name__ == "__main__":
    main()
