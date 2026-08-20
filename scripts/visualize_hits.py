#!/usr/bin/env python3
"""Render a per-gene PNG+SVG showing all kmerseek hits mapped onto the query protein.

For one query sequence, draws its full length as a bar, overlays every matched
target's actual matched regions positioned to scale (stacked into lanes when
hits overlap), and prints the query / encoded-alphabet / target alignment
beneath each hit. Hits are numbered instead of connected to their alignment
block with leader lines, since lines cross when hits interleave.

Input is the CSV produced by `kmerseek search -o results.csv` (one row per
matched region: region_start, region_end, region_subseq, target_start, target_end,
target_subseq, moltype_seq, ...) plus the query FASTA, used to draw the full-length
protein bar and to get the exact query names to plot.

Usage:
    python visualize_hits.py --csv results.csv --query-fasta query.fasta \
        --output-dir hits_png/
    python visualize_hits.py --csv results.csv --query-fasta query.fasta \
        --output-dir hits_png/ --query-name CED9
"""

import argparse
import os
import re
import textwrap
from collections import defaultdict

import matplotlib
import polars as pl

matplotlib.use("Agg")
# Keep SVG text as real <text> elements (selectable/copyable), not vector outlines.
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Categorical colors always come from a built-in matplotlib qualitative colormap,
# never a hand-picked palette. Picked by how many distinct targets need distinct
# colors: small counts get the punchier 8-color sets, larger counts step up to
# maps with more distinguishable steps (repeating past 20, the largest available).
_QUALITATIVE_CMAPS = [(8, "Set2"), (10, "tab10"), (12, "Set3"), (20, "tab20")]


def get_palette(n_colors):
    for max_n, cmap_name in _QUALITATIVE_CMAPS:
        if n_colors <= max_n:
            return list(matplotlib.colormaps[cmap_name].colors)
    return list(matplotlib.colormaps["tab20"].colors)

INK = "#0b0b0b"
SECONDARY_INK = "#52514e"
MUTED = "#898781"
BASELINE = "#c3c2b7"
BAR_FILL = "#e1e0d9"
SURFACE = "#fcfcfb"

SEQ_WRAP_WIDTH = 60
MAX_SEQ_LINES = 3


def _relative_luminance(rgb):
    """Perceptual luminance of an (r, g, b[, a]) tuple in 0..1, for choosing
    readable (white vs. dark) text on top of an arbitrary palette color."""
    r, g, b = rgb[:3]
    return 0.299 * r + 0.587 * g + 0.114 * b


def read_fasta_lengths(fasta_path):
    """Return {full_header: sequence_length}, keyed the same way needletail's
    record.id() names a query (the whole header line after '>')."""
    lengths = {}
    name = None
    length = 0
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:]
                length = 0
            else:
                length += len(line.strip())
        if name is not None:
            lengths[name] = length
    return lengths


def short_label(name, max_len=28):
    """Shorten a pipe-delimited header to a display id; falls back to the
    first word. Two header shapes are recognized:
    - GENCODE ('ENSP...|ENST...|ENSG...|OTTHUMG...|OTTHUMT...|GENE-201|GENE|len'):
      the gene symbol is field index 6, not 2 -- using index 2 would show the
      ENSG accession instead of a readable gene name.
    - UniProt-style ('sp|P10415|BCL2_HUMAN Apoptosis...'): the name is field
      index 2.
    """
    parts = name.split("|")
    if len(parts) >= 8:
        short = parts[6].split(" ")[0]
    elif len(parts) >= 3:
        short = parts[2].split(" ")[0]
    else:
        short = name.split(" ")[0]
    if len(short) > max_len:
        short = short[: max_len - 1] + "…"
    return short


def safe_filename(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", short_label(name, max_len=60)).strip("_")


def resolve_query_names(query_name_arg, all_query_names):
    """Match --query-name against the full FASTA headers found in the CSV.

    Accepts the exact header (backward compatible), but also a short,
    case-insensitive substring like "CED9" or "ced9_caeel" matched against
    either the full header or its short_label() -- typing the whole
    'sp|P41958|CED9_CAEEL Apoptosis regulator...' header is not required.
    """
    if query_name_arg is None:
        return sorted(all_query_names)
    if query_name_arg in all_query_names:
        return [query_name_arg]
    needle = query_name_arg.lower()
    matches = [
        name for name in all_query_names
        if needle in name.lower() or needle in short_label(name, max_len=len(name)).lower()
    ]
    return sorted(set(matches))


def scan_csv(csv_path):
    """Open the results CSV lazily -- nothing is read off disk until a query
    is `.collect()`-ed, so the whole file never has to fit in memory at once.
    infer_schema_length=None scans every row for dtypes up front, since a
    p-value column of mostly-0 values with rare 1e-300-style outliers
    further down the file could otherwise get mis-inferred as a narrower type."""
    return pl.scan_csv(csv_path, infer_schema_length=None)


def load_query_names(lazy_rows):
    """Every distinct query_name in the CSV. Only the query_name column is
    read to compute this, not the whole (much wider) row."""
    return set(lazy_rows.select("query_name").unique().collect().get_column("query_name"))


def load_rows_for_queries(lazy_rows, query_names):
    """Collect only the rows belonging to query_names, in one pass. Filtering
    happens before collecting, so the common `--query-name` case (a small
    subset) never materializes the rest of the CSV."""
    return lazy_rows.filter(pl.col("query_name").is_in(list(query_names))).collect()


def iter_query_rows(df):
    """Yield (query_name, rows) one query at a time -- rows as plain dicts,
    matching the shape the rest of this module already expects from
    csv.DictReader. A single grouped pass over the collected DataFrame, so
    only one query's rows are ever converted to Python dicts at a time,
    instead of the whole result up front (partition_by(as_dict=True) does
    that eagerly, roughly doubling peak memory for a "plot every query" run).
    Re-filtering the DataFrame once per query instead of grouping would avoid
    that too, but stays O(queries x rows) -- no better than the old
    pure-Python approach once more than a handful of queries are requested."""
    for (query_name,), group in df.group_by("query_name", maintain_order=True):
        yield query_name, group.to_dicts()


def merge_regions_by_target(rows, gap_merge):
    """Group a query's matched regions by target, merging regions that are
    within `gap_merge` residues of each other (in query coordinates) into a
    single hit. Returns one dict per hit, keeping every individual region's
    real position -- a hit's `start`/`end` is only the outer envelope, not a
    claim that the whole span matched."""
    by_target = defaultdict(list)
    for row in rows:
        by_target[row["target_name"]].append(row)

    hits = []
    for target_name, target_rows in by_target.items():
        target_rows.sort(key=lambda r: int(r["region_start"]))
        cluster = [target_rows[0]]
        for row in target_rows[1:]:
            if int(row["region_start"]) - int(cluster[-1]["region_end"]) <= gap_merge:
                cluster.append(row)
            else:
                hits.append(_build_hit(target_name, cluster))
                cluster = [row]
        hits.append(_build_hit(target_name, cluster))
    return hits


def _union_coverage(regions):
    """Total residues covered by the union of (start, end) intervals -- summing
    raw lengths would double-count overlapping sliding-window matches."""
    merged_end = None
    covered = 0
    for start, end in sorted(regions):
        if merged_end is None or start > merged_end:
            covered += end - start
            merged_end = end
        elif end > merged_end:
            covered += end - merged_end
            merged_end = end
    return covered


def _target_best_rows(rows):
    """{target_name: the CSV row of its best (highest-scoring) region}, one entry per
    distinct target.

    region_poisson_score is -log10 of the region's Poisson tail probability, so
    bigger means more surprising (see MatchedRegion::poisson_score in
    src/rust/search.rs). It's a ranking heuristic, not a calibrated p-value --
    but it is still the number that decides which hits this plot surfaces, so
    it is what gets corrected here, not the whole-query p-value. This plot
    draws regions, and `kmerseek search` now reports a hit when either scope
    clears, so a real sub-protein domain call routinely carries an unimpressive
    whole-query p-value (BCL2/CED9: 0.99 whole-query, score ~3.16 i.e. p=0.0007
    for its region). Correcting the whole-query number instead would push
    those hits to q~1 and let --max-hits cut them, hiding what the region
    scoring exists to surface.

    Unlike the query-level stats, region scores differ row to row, so take the
    strongest (highest-scoring) region as the target's evidence."""
    best = {}
    for row in rows:
        name = row["target_name"]
        if name not in best or float(row["region_poisson_score"]) > float(best[name]["region_poisson_score"]):
            best[name] = row
    return best


def _target_scores(rows):
    """{target_name: best (highest) region score}, one entry per distinct target."""
    return {name: float(row["region_poisson_score"]) for name, row in _target_best_rows(rows).items()}


def _target_tail_probabilities(rows):
    """{target_name: raw Poisson tail probability of its best-scoring region}, one entry per
    distinct target.

    Reads region_tail_probability straight off the same row _target_scores takes its score
    from (both columns come from the same MatchedRegion), instead of reconstructing a
    probability by undoing the -log10 transform on the score. benjamini_hochberg needs real
    probabilities to correct, not scores."""
    return {name: float(row["region_tail_probability"]) for name, row in _target_best_rows(rows).items()}


def benjamini_hochberg(pvalues):
    """Benjamini-Hochberg FDR-corrected p-values (q-values) for multiple-testing
    correction across every target tested for one query. pvalues: {key: p}.
    Returns {key: q}, monotonic non-decreasing as p increases, capped at 1.0."""
    m = len(pvalues)
    ranked = sorted(pvalues.items(), key=lambda kv: kv[1])
    corrected = {}
    min_so_far = 1.0
    for rank in range(m, 0, -1):
        key, p = ranked[rank - 1]
        min_so_far = min(min_so_far, p * m / rank, 1.0)
        corrected[key] = min_so_far
    return corrected


def _build_hit(target_name, cluster):
    region_rows = sorted(cluster, key=lambda r: int(r["region_start"]))
    regions = [(int(r["region_start"]), int(r["region_end"])) for r in region_rows]
    return {
        "target_name": target_name,
        "start": min(r[0] for r in regions),
        "end": max(r[1] for r in regions),
        "regions": regions,
        "region_rows": region_rows,
        "coverage": _union_coverage(regions),
        "n_regions": len(cluster),
        # containment/jaccard/query_* are query-target *result* stats, not
        # per-region -- every row in the cluster carries the same value (it's the
        # same query-target pair); max() is just a safe pick.
        "containment": max(float(r["containment"]) for r in cluster),
        "jaccard": max(float(r["jaccard"]) for r in cluster),
        "query_enrichment": max(float(r["query_enrichment"]) for r in cluster),
        "query_poisson_pvalue": max(float(r["query_poisson_pvalue"]) for r in cluster),
        # Region stats do vary per row; the strongest (highest-scoring) region is the hit's
        # evidence.
        "region_poisson_score": max(float(r["region_poisson_score"]) for r in cluster),
        "region_enrichment": max(float(r["region_enrichment"]) for r in cluster),
        "moltype": region_rows[0]["moltype"],
    }


def assign_lanes(hits):
    """Greedy interval scheduling: put each hit in the first lane whose last
    hit doesn't overlap it, so overlapping hits stack into separate rows."""
    lanes = []  # lanes[i] = end position of the last hit placed in lane i
    for hit in sorted(hits, key=lambda h: h["start"]):
        placed = False
        for lane_idx, lane_end in enumerate(lanes):
            if hit["start"] >= lane_end:
                hit["lane"] = lane_idx
                lanes[lane_idx] = hit["end"]
                placed = True
                break
        if not placed:
            hit["lane"] = len(lanes)
            lanes.append(hit["end"])
    return len(lanes)


def wrap_seq(seq):
    lines = textwrap.wrap(seq, SEQ_WRAP_WIDTH) or [""]
    if len(lines) > MAX_SEQ_LINES:
        lines = lines[:MAX_SEQ_LINES]
        lines[-1] = lines[-1][: SEQ_WRAP_WIDTH - 1] + "…"
    return lines


class GenePlot:
    """Renders one gene's PNG+SVG: a full-length protein bar, every hit's real
    matched regions positioned to scale below it (numbered, stacked into lanes
    when hits overlap), and an alignment block per hit underneath.

    Usage: GenePlot(query_name, query_length, hits).render(output_paths).
    """

    # --- layout constants (inches unless noted) ---
    TRACK_H = 0.35  # height of one hit lane / the protein bar
    LANE_GAP = 0.08  # vertical gap between stacked lanes, and between the bar and lane 0
    MIN_REGION_FRACTION = 0.004  # a region always draws at least this fraction of query_length wide
    FIG_W = 11
    BADGE_FONTSIZE = 7.5
    LINE_H = 0.16  # height of one line of sequence/label text in the alignment blocks
    LABEL_X_OFFSET = 0.022  # left indent of a header/region-label/group-label line
    SEQ_X_OFFSET = 0.09  # left indent of sequence text, past the "query:"/"hp:"/"target:" label
    HEADER_H = 0.46  # space for a hit's name/span line + stats line + padding before its first region
    GROUP_GAP = 0.20  # extra gap (in LINE_H units) after each query/moltype/target group
    HIT_GAP = 0.15  # gap after a hit's last region, before the next hit's header

    def __init__(self, query_name, query_length, hits, dpi=200):
        self.query_name = query_name
        self.query_length = query_length
        self.dpi = dpi
        self.n_lanes = assign_lanes(hits)
        self.hits = sorted(hits, key=lambda h: h["start"])
        self.colors = self._assign_hit_colors()
        self.text_height_in = max(sum(self._hit_text_height(h) for h in self.hits), 0.1)
        self.fig = self.ax = self.ax_text = None

    def render(self, output_paths):
        """Draw the full figure and save it to every path in output_paths."""
        self._make_figure()
        self._draw_track()
        self._draw_alignment_blocks()
        for output_path in output_paths:
            self.fig.savefig(output_path, dpi=self.dpi, facecolor=SURFACE, bbox_inches="tight")
        plt.close(self.fig)
        return self

    def _assign_hit_colors(self):
        """One color per target, not per hit-fragment -- a target split into
        multiple hits (e.g. by a gap too large to merge) keeps one color
        throughout. Tie-break on name: iterating a set is hash-seed dependent,
        so without it, two targets tied on start position could swap colors
        between runs. Palette size scales with how many targets need distinguishing."""
        target_order = sorted(
            {h["target_name"] for h in self.hits},
            key=lambda t: (min(h["start"] for h in self.hits if h["target_name"] == t), t),
        )
        palette = get_palette(len(target_order))
        target_colors = {t: palette[i % len(palette)] for i, t in enumerate(target_order)}
        return {id(hit): target_colors[hit["target_name"]] for hit in self.hits}

    def _hit_text_height(self, hit):
        """Vertical space (inches) a hit's alignment block will take, including
        every one of its regions -- not just a single representative one."""
        height = self.HEADER_H
        show_region_labels = hit["n_regions"] > 1
        for row in hit["region_rows"]:
            if show_region_labels:
                height += self.LINE_H
            for key in ("region_subseq", "moltype_seq", "target_subseq"):
                height += len(wrap_seq(row[key])) * self.LINE_H + self.LINE_H * self.GROUP_GAP
        return height + self.HIT_GAP

    def _make_figure(self):
        track_height_in = 0.9 + self.n_lanes * (self.TRACK_H + self.LANE_GAP)
        fig_h = track_height_in + self.text_height_in + 0.6
        self.fig, (self.ax, self.ax_text) = plt.subplots(
            2, 1, figsize=(self.FIG_W, fig_h), dpi=self.dpi,
            gridspec_kw={"height_ratios": [track_height_in, self.text_height_in], "hspace": 0.35},
        )
        self.fig.patch.set_facecolor(SURFACE)
        self.ax.set_facecolor(SURFACE)
        self.ax_text.set_facecolor(SURFACE)

    # --- track (protein bar + hit lanes + number badges) ---

    def _draw_track(self):
        bar_y = self.n_lanes * (self.TRACK_H + self.LANE_GAP)
        self._draw_protein_bar(bar_y)
        self._style_track_axes(bar_y)
        widest_region = self._draw_hit_regions()
        self._draw_hit_numbers(widest_region)

    def _draw_protein_bar(self, bar_y):
        ax = self.ax
        ax.add_patch(Rectangle((0, bar_y), self.query_length, self.TRACK_H,
                                facecolor=BAR_FILL, edgecolor=BASELINE, linewidth=1))
        ax.text(0, bar_y + self.TRACK_H + 0.08, "1", ha="left", va="bottom", fontsize=8, color=MUTED)
        ax.text(self.query_length, bar_y + self.TRACK_H + 0.08, f"{self.query_length}aa",
                ha="right", va="bottom", fontsize=8, color=MUTED)
        ax.text(0, bar_y + self.TRACK_H + 0.28, short_label(self.query_name, max_len=60),
                ha="left", va="bottom", fontsize=12, color=INK, fontweight="bold")

    def _style_track_axes(self, bar_y):
        ax = self.ax
        ax.set_xlim(-self.query_length * 0.02, self.query_length * 1.02)
        ax.set_ylim(-0.1, bar_y + self.TRACK_H + 0.6)
        ax.set_yticks([])
        ax.set_xlabel("Query position (aa)", fontsize=9, color=SECONDARY_INK)
        ax.tick_params(axis="x", colors=MUTED, labelsize=8)
        for spine in ("top", "left", "right"):
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color(BASELINE)

    def _draw_hit_regions(self):
        """Draw every hit's real matched regions, linked by a thin envelope line
        across their span -- never implying more coverage than there is. Returns
        each hit's widest drawn box (for number placement), keyed by id(hit)."""
        min_region_w = self.query_length * self.MIN_REGION_FRACTION
        widest_region = {}
        for hit in self.hits:
            color = self.colors[id(hit)]
            y = hit["lane"] * (self.TRACK_H + self.LANE_GAP)
            if hit["end"] - hit["start"] > min_region_w:
                self.ax.plot([hit["start"], hit["end"]], [y + self.TRACK_H / 2, y + self.TRACK_H / 2],
                             color=color, linewidth=1.2, alpha=0.55, zorder=1, solid_capstyle="butt")
            widest_region[id(hit)] = self._draw_hit_boxes(hit, y, color, min_region_w)
        return widest_region

    def _draw_hit_boxes(self, hit, y, color, min_region_w):
        """Draw every one of a hit's regions as a solid box; return the (start,
        width) of its widest drawn box, used to place that hit's number badge."""
        widest = None
        for r_start, r_end in hit["regions"]:
            width = max(r_end - r_start, min_region_w)
            self.ax.add_patch(Rectangle((r_start, y), width, self.TRACK_H,
                                         facecolor=color, edgecolor=SURFACE, linewidth=0.6, zorder=2))
            if widest is None or width > widest[1]:
                widest = (r_start, width)
        return widest

    def _draw_hit_numbers(self, widest_region):
        """Number badges go *inside* each hit's widest box, never floating above
        it where tight lane spacing could make them collide with the lane above."""
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()
        for number, hit in enumerate(self.hits, start=1):
            y = hit["lane"] * (self.TRACK_H + self.LANE_GAP)
            r_start, width = widest_region[id(hit)]
            self._draw_one_number_badge(renderer, number, r_start, y, width, self.colors[id(hit)])

    def _draw_one_number_badge(self, renderer, number, r_start, y, width, color):
        """Draw a hit's number centered in its box; remove it if it doesn't fit --
        a skipped number (still shown in the alignment block below) beats an
        overflowing one that collides with whatever is next to it."""
        text_color = "white" if _relative_luminance(color) < 0.55 else INK
        label = self.ax.text(r_start + width / 2, y + self.TRACK_H / 2, str(number),
                              ha="center", va="center", fontsize=self.BADGE_FONTSIZE,
                              color=text_color, fontweight="bold", zorder=3)
        bbox = label.get_window_extent(renderer=renderer)
        (x0, y0), (x1, y1) = self.ax.transData.transform(
            [(r_start, y), (r_start + width, y + self.TRACK_H)])
        if bbox.width > x1 - x0 - 2 or bbox.height > y1 - y0 - 2:
            label.remove()

    # --- alignment blocks (header + every region's query/hp/target lines) ---

    def _draw_alignment_blocks(self):
        """Lay out every hit's alignment block top-to-bottom. ax_text has its own
        unit-free grid so text offsets don't get squashed by the aa-position scale
        used in the track axes; each block's height is exactly _hit_text_height."""
        self.ax_text.set_xlim(0, 1)
        self.ax_text.set_ylim(-self.text_height_in, 0.15)
        self.ax_text.axis("off")
        y_cursor = 0.0
        for i, hit in enumerate(self.hits):
            self._draw_hit_alignment_block(i + 1, hit, self.colors[id(hit)], y_cursor)
            y_cursor -= self._hit_text_height(hit)

    def _draw_hit_alignment_block(self, index, hit, color, y_top):
        """Draw one hit's header (name/span line + a stats line), then every
        region's alignment (not just one representative) -- e.g. a 2-region hit
        like BAK_HUMAN shows both."""
        self.ax_text.text(self.LABEL_X_OFFSET, y_top, self._hit_header_text(index, hit),
                           ha="left", va="top", fontsize=8, color=color,
                           fontweight="bold", family="sans-serif")
        self.ax_text.text(self.LABEL_X_OFFSET, y_top - self.LINE_H, self._hit_stats_text(hit),
                           ha="left", va="top", fontsize=7, color=MUTED, family="sans-serif")
        y = y_top - self.HEADER_H
        show_region_labels = hit["n_regions"] > 1
        for r_idx, row in enumerate(hit["region_rows"], start=1):
            if show_region_labels:
                self.ax_text.text(self.LABEL_X_OFFSET, y, self._region_label_text(r_idx, hit, row),
                                   ha="left", va="top", fontsize=7, color=MUTED, family="sans-serif")
                y -= self.LINE_H
            y = self._draw_region_alignment(row, y)

    @staticmethod
    def _hit_header_text(index, hit):
        span = hit["end"] - hit["start"]
        extra = f", {hit['n_regions']} regions" if hit["n_regions"] > 1 else ""
        return (f"{index}. {short_label(hit['target_name'])}  "
                f"(covers {hit['coverage']}/{span}aa span at {hit['start'] + 1}-{hit['end']}aa{extra})")

    @staticmethod
    def _hit_stats_text(hit):
        # corrected_pvalue is attached by _render_query (a BH-FDR q-value across
        # every target tested for this query); absent when a hit is built and
        # plotted directly, e.g. in tests, without going through that pipeline.
        q_value = hit.get("corrected_pvalue")
        q_part = f"   region q={q_value:.2g}" if q_value is not None else ""
        # Both scopes are shown: either one can be what got this hit reported, and
        # seeing them side by side is how you tell a whole-protein match from a
        # localized domain call.
        return (f"containment={hit['containment']:.2f}   jaccard={hit['jaccard']:.3f}   "
                f"region enrich={hit['region_enrichment']:.2f}   "
                f"region score={hit['region_poisson_score']:.2g}{q_part}   "
                f"query p={hit['query_poisson_pvalue']:.2g}")

    @staticmethod
    def _region_label_text(r_idx, hit, row):
        r_start, r_end = int(row["region_start"]), int(row["region_end"])
        return (f"region {r_idx}/{hit['n_regions']}:  {r_start + 1}-{r_end}aa, "
                f"score={float(row['region_poisson_score']):.2g}")

    def _draw_region_alignment(self, row, y):
        """Draw one region's query/moltype/target lines; return the y cursor after it."""
        groups = (("query", row["region_subseq"]), (row["moltype"], row["moltype_seq"]),
                  ("target", row["target_subseq"]))
        for label, seq in groups:
            self.ax_text.text(self.LABEL_X_OFFSET, y, f"{label}:", ha="left", va="top",
                               fontsize=7.5, color=MUTED, family="monospace")
            for line in wrap_seq(seq):
                self.ax_text.text(self.SEQ_X_OFFSET, y, line, ha="left", va="top",
                                   fontsize=7.5, color=INK, family="monospace")
                y -= self.LINE_H
            y -= self.LINE_H * self.GROUP_GAP
        return y


def plot_gene(query_name, query_length, hits, output_paths, dpi=200):
    """Render one gene's PNG+SVG. Thin functional wrapper around GenePlot."""
    GenePlot(query_name, query_length, hits, dpi=dpi).render(output_paths)


def _build_arg_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", required=True, help="kmerseek search results CSV")
    parser.add_argument("--query-fasta", required=True, help="Query FASTA (for full-length protein bars)")
    parser.add_argument("--query-name", default=None,
                         help="Gene to plot -- a short case-insensitive substring like 'CED9' matched against "
                              "the FASTA header (the exact full header also works). Default: plot every query "
                              "found in the CSV.")
    parser.add_argument("--output-dir", required=True, help="Directory to write one PNG+SVG pair per gene into")
    parser.add_argument("--gap-merge", type=int, default=10,
                         help="Merge same-target regions within this many residues into one hit (default: 10)")
    parser.add_argument("--min-containment", type=float, default=0.0,
                         help="Drop matched regions below this containment before plotting (default: 0.0, keep all)")
    parser.add_argument("--max-hits", type=int, default=None,
                         help="Keep only the top N distinct targets per gene, ranked by each target's "
                              "BH-corrected q-value, most significant first (default: unlimited); all of a kept "
                              "target's hit spans are still shown. Proteome-scale searches can produce dozens of "
                              "targets per gene, making the default unlimited figure very tall -- use this to cap it.")
    parser.add_argument("--dpi", type=int, default=200)
    return parser


def _cap_to_top_targets(hits, max_hits, corrected_pvalues):
    """Keep only hits belonging to the top max_hits distinct targets, ranked by
    each target's BH-corrected q-value (most significant first, ties broken by
    containment) -- capping by hit-span count instead would let one
    heavily-fragmented target crowd out every other target."""
    target_names = {h["target_name"] for h in hits}
    if max_hits is None or len(target_names) <= max_hits:
        return hits
    best_containment = {}
    for h in hits:
        best_containment[h["target_name"]] = max(
            best_containment.get(h["target_name"], 0.0), h["containment"])
    top_targets = set(sorted(
        target_names,
        key=lambda t: (corrected_pvalues[t], -best_containment[t]),
    )[:max_hits])
    return [h for h in hits if h["target_name"] in top_targets]


def _render_query(query_name, query_rows, query_length, args):
    """Build one gene's hits from its CSV rows and render its PNG+SVG pair.

    query_rows is every row for this query in the CSV, unfiltered by
    --min-containment -- the BH-FDR q-value must be corrected across every
    target actually tested, not just the ones that end up displayed, or it
    would understate how many comparisons were made.
    """
    corrected_pvalues = benjamini_hochberg(_target_tail_probabilities(query_rows))
    display_rows = [r for r in query_rows if float(r["containment"]) >= args.min_containment]
    if not display_rows:
        print(f"Skipping '{query_name}': no hits above --min-containment {args.min_containment}")
        return

    hits = _cap_to_top_targets(
        merge_regions_by_target(display_rows, args.gap_merge), args.max_hits, corrected_pvalues)
    for hit in hits:
        hit["corrected_pvalue"] = corrected_pvalues[hit["target_name"]]

    base = os.path.join(args.output_dir, f"{safe_filename(query_name)}.hits")
    plot_gene(query_name, query_length, hits, [f"{base}.png", f"{base}.svg"], dpi=args.dpi)
    n_targets = len({h["target_name"] for h in hits})
    print(f"Wrote {base}.png / .svg ({n_targets} targets, {len(hits)} hits, {len(display_rows)} regions)")


def main():
    args = _build_arg_parser().parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    lengths = read_fasta_lengths(args.query_fasta)
    lazy_rows = scan_csv(args.csv)

    all_query_names = load_query_names(lazy_rows)
    query_names = resolve_query_names(args.query_name, all_query_names)
    if args.query_name is not None and not query_names:
        available = ", ".join(sorted(short_label(n) for n in all_query_names))
        print(f"No query matching '{args.query_name}' found in {args.csv}. Available: {available}")
        return

    df = load_rows_for_queries(lazy_rows, query_names)
    for query_name, query_rows in iter_query_rows(df):
        if query_name not in lengths:
            print(f"Skipping '{query_name}': not found in {args.query_fasta}")
            continue
        _render_query(query_name, query_rows, lengths[query_name], args)


if __name__ == "__main__":
    main()
