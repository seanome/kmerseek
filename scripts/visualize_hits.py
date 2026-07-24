#!/usr/bin/env python3
"""Render a per-gene PNG+SVG showing all kmerseek hits mapped onto the query protein.

For one query sequence, draws its full length as a bar, overlays every matched
target's actual matched regions positioned to scale (stacked into lanes when
hits overlap), and prints the query / encoded-alphabet / target alignment
beneath each hit. Hits are numbered instead of connected to their alignment
block with leader lines, since lines cross when hits interleave.

Input is the CSV produced by `kmerseek search -o results.csv` (one row per
matched region: query_start, query_end, query_subseq, target_start, target_end,
target_subseq, moltype_seq, ...) plus the query FASTA, used to draw the full-length
protein bar and to get the exact query names to plot.

Usage:
    python visualize_hits.py --csv results.csv --query-fasta query.fasta \
        --output-dir hits_png/
    python visualize_hits.py --csv results.csv --query-fasta query.fasta \
        --output-dir hits_png/ --query-name CED9
"""

import argparse
import csv
import os
import re
import textwrap
from collections import defaultdict

import matplotlib

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
    """Shorten a UniProt-style header ('sp|P10415|BCL2_HUMAN Apoptosis...')
    to a display id ('BCL2_HUMAN'); falls back to the first word."""
    parts = name.split("|")
    short = parts[2].split(" ")[0] if len(parts) >= 3 else name.split(" ")[0]
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


def load_rows(csv_path):
    with open(csv_path, newline="") as fh:
        return list(csv.DictReader(fh))


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
        target_rows.sort(key=lambda r: int(r["query_start"]))
        cluster = [target_rows[0]]
        for row in target_rows[1:]:
            if int(row["query_start"]) - int(cluster[-1]["query_end"]) <= gap_merge:
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


def _build_hit(target_name, cluster):
    region_rows = sorted(cluster, key=lambda r: int(r["query_start"]))
    regions = [(int(r["query_start"]), int(r["query_end"])) for r in region_rows]
    return {
        "target_name": target_name,
        "start": min(r[0] for r in regions),
        "end": max(r[1] for r in regions),
        "regions": regions,
        "region_rows": region_rows,
        "coverage": _union_coverage(regions),
        "n_regions": len(cluster),
        "containment": max(float(r["containment"]) for r in cluster),
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


HEADER_H = 0.30  # space for a hit's header line + padding before its first region
GROUP_GAP = 0.20  # extra gap (in line_h units) after each query/moltype/target group
HIT_GAP = 0.15  # gap after a hit's last region, before the next hit's header


def _hit_text_height(hit, line_h):
    """Total vertical space (inches) a hit's alignment block will take, including
    every one of its regions -- not just a single representative one."""
    height = HEADER_H
    show_region_labels = hit["n_regions"] > 1
    for row in hit["region_rows"]:
        if show_region_labels:
            height += line_h
        for key in ("query_subseq", "moltype_seq", "target_subseq"):
            height += len(wrap_seq(row[key])) * line_h + line_h * GROUP_GAP
    return height + HIT_GAP


def plot_gene(query_name, query_length, hits, output_paths, dpi=200):
    n_lanes = assign_lanes(hits)
    hits_by_start = sorted(hits, key=lambda h: h["start"])

    # Color follows the target, not the fragment: a target split into multiple
    # hits (e.g. by a gap too large to merge) keeps one color throughout. Palette
    # size scales with how many distinct targets need distinguishing.
    target_order = sorted({h["target_name"] for h in hits_by_start},
                           key=lambda t: min(h["start"] for h in hits_by_start if h["target_name"] == t))
    palette = get_palette(len(target_order))
    target_colors = {t: palette[i % len(palette)] for i, t in enumerate(target_order)}
    colors = {id(hit): target_colors[hit["target_name"]] for hit in hits_by_start}

    # One alignment block per hit, one row per line -- a 2-column grid let long
    # headers (variable-length target names + stats) spill into the next column.
    # Every region in a hit gets its own alignment shown, not just one
    # representative, so a block's height depends on how many regions it has.
    line_h = 0.16
    text_height_in = max(sum(_hit_text_height(h, line_h) for h in hits_by_start), 0.1)

    track_h = 0.35
    lane_gap = 0.08
    track_height_in = 0.9 + n_lanes * (track_h + lane_gap)
    fig_h = track_height_in + text_height_in + 0.6
    fig_w = 11

    fig, (ax, ax_text) = plt.subplots(
        2, 1, figsize=(fig_w, fig_h), dpi=dpi,
        gridspec_kw={"height_ratios": [track_height_in, text_height_in], "hspace": 0.35},
    )
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    ax_text.set_facecolor(SURFACE)

    # --- full-length protein bar ---
    # Same (track_h + lane_gap) stacking step used between hit lanes, so the gap
    # above the topmost hit lane doesn't end up wider than the gaps between lanes.
    bar_y = n_lanes * (track_h + lane_gap)
    ax.add_patch(
        Rectangle(
            (0, bar_y), query_length, track_h,
            facecolor=BAR_FILL, edgecolor=BASELINE, linewidth=1,
        )
    )
    ax.text(0, bar_y + track_h + 0.08, "1", ha="left", va="bottom",
             fontsize=8, color=MUTED)
    ax.text(query_length, bar_y + track_h + 0.08, f"{query_length}aa",
             ha="right", va="bottom", fontsize=8, color=MUTED)
    ax.text(0, bar_y + track_h + 0.28, short_label(query_name, max_len=60),
             ha="left", va="bottom", fontsize=12, color=INK, fontweight="bold")

    ax.set_xlim(-query_length * 0.02, query_length * 1.02)
    ax.set_ylim(-0.1, bar_y + track_h + 0.6)
    ax.set_yticks([])
    ax.set_xlabel("Query position (aa)", fontsize=9, color=SECONDARY_INK)
    ax.tick_params(axis="x", colors=MUTED, labelsize=8)
    for spine in ("top", "left", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(BASELINE)

    # --- hit lanes: a thin envelope line spans the hit's full footprint (so
    # nearby fragments read as one hit), with the *actual* matched regions
    # drawn as solid blocks on top -- never implying more coverage than there is.
    min_region_w = query_length * 0.004
    widest_region = {}  # id(hit) -> (r_start, drawn_width) of its widest drawn box
    for hit in hits_by_start:
        color = colors[id(hit)]
        y = hit["lane"] * (track_h + lane_gap)

        if hit["end"] - hit["start"] > min_region_w:
            ax.plot([hit["start"], hit["end"]], [y + track_h / 2, y + track_h / 2],
                     color=color, linewidth=1.2, alpha=0.55, zorder=1, solid_capstyle="butt")

        for r_start, r_end in hit["regions"]:
            width = max(r_end - r_start, min_region_w)
            ax.add_patch(
                Rectangle(
                    (r_start, y), width, track_h,
                    facecolor=color, edgecolor=SURFACE, linewidth=0.6, zorder=2,
                )
            )
            best = widest_region.get(id(hit))
            if best is None or width > best[1]:
                widest_region[id(hit)] = (r_start, width)

    # Numbers go *inside* each hit's widest box (never floating above it, where
    # tight lane spacing let them collide with the box in the lane above). Drawn
    # after every box so we can measure each number's actual rendered footprint
    # and skip it outright -- rather than let it overflow -- when the box is too
    # small; the number is still available via the alignment block below.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for number, hit in enumerate(hits_by_start, start=1):
        color = colors[id(hit)]
        y = hit["lane"] * (track_h + lane_gap)
        r_start, width = widest_region[id(hit)]
        text_color = "white" if _relative_luminance(color) < 0.55 else INK
        label = ax.text(r_start + width / 2, y + track_h / 2, str(number),
                         ha="center", va="center", fontsize=7.5, color=text_color,
                         fontweight="bold", zorder=3)
        bbox = label.get_window_extent(renderer=renderer)
        box_bbox = ax.transData.transform([(r_start, y), (r_start + width, y + track_h)])
        box_width_px = box_bbox[1][0] - box_bbox[0][0]
        box_height_px = box_bbox[1][1] - box_bbox[0][1]
        if bbox.width > box_width_px - 2 or bbox.height > box_height_px - 2:
            label.remove()

    # --- alignment text blocks, in reading order top-to-bottom, one per hit ---
    # Every region belonging to a hit is shown (not just one representative),
    # so a multi-region hit (e.g. BAK_HUMAN with 2 separate matches) shows both.
    # ax_text has its own unit-free grid (0..1 wide) so text offsets don't get
    # squashed by the aa-position scale used in `ax`; height is the exact sum
    # of each hit's content, computed once above via _hit_text_height.
    ax_text.set_xlim(0, 1)
    ax_text.set_ylim(-text_height_in, 0.15)
    ax_text.axis("off")

    label_x_offset = 0.022
    seq_x_offset = 0.09
    y_cursor = 0.0
    for i, hit in enumerate(hits_by_start):
        x0 = 0
        y_top = y_cursor
        color = colors[id(hit)]
        n_regions = hit["n_regions"]
        show_region_labels = n_regions > 1

        header = (
            f"{i + 1}. {short_label(hit['target_name'])}  "
            f"(covers {hit['coverage']}/{hit['end'] - hit['start']}aa span "
            f"at {hit['start']+1}-{hit['end']}aa, "
            f"containment={hit['containment']:.2f}"
            + (f", {n_regions} regions" if n_regions > 1 else "")
            + ")"
        )
        ax_text.text(x0 + label_x_offset, y_top, header, ha="left", va="top",
                      fontsize=8, color=color, fontweight="bold", family="sans-serif")

        y = y_top - HEADER_H
        for r_idx, row in enumerate(hit["region_rows"], start=1):
            if show_region_labels:
                r_start, r_end = int(row["query_start"]), int(row["query_end"])
                region_label = (
                    f"region {r_idx}/{n_regions}:  {r_start + 1}-{r_end}aa, "
                    f"containment={float(row['containment']):.2f}"
                )
                ax_text.text(x0 + label_x_offset, y, region_label, ha="left", va="top",
                              fontsize=7, color=MUTED, family="sans-serif")
                y -= line_h

            q_lines = wrap_seq(row["query_subseq"])
            m_lines = wrap_seq(row["moltype_seq"])
            t_lines = wrap_seq(row["target_subseq"])
            for label, lines in (("query", q_lines), (row["moltype"], m_lines), ("target", t_lines)):
                ax_text.text(x0 + label_x_offset, y, f"{label}:", ha="left", va="top",
                              fontsize=7.5, color=MUTED, family="monospace")
                for line in lines:
                    ax_text.text(x0 + seq_x_offset, y, line, ha="left", va="top",
                                  fontsize=7.5, color=INK, family="monospace")
                    y -= line_h
                y -= line_h * GROUP_GAP

        y_cursor = y_top - _hit_text_height(hit, line_h)

    for output_path in output_paths:
        fig.savefig(output_path, dpi=dpi, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


def main():
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
                         help="Keep only the top N distinct targets per gene, ranked by each target's best "
                              "containment (default: unlimited); all of a kept target's hit spans are still shown. "
                              "Proteome-scale searches can produce dozens of targets per gene, making the default "
                              "unlimited figure very tall -- use this to cap it.")
    parser.add_argument("--dpi", type=int, default=200)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    lengths = read_fasta_lengths(args.query_fasta)
    rows = load_rows(args.csv)
    if args.min_containment > 0.0:
        rows = [r for r in rows if float(r["containment"]) >= args.min_containment]

    all_query_names = {r["query_name"] for r in rows}
    query_names = resolve_query_names(args.query_name, all_query_names)
    if args.query_name is not None and not query_names:
        available = ", ".join(sorted(short_label(n) for n in all_query_names))
        print(f"No query matching '{args.query_name}' found in {args.csv}. "
              f"Available: {available}")
        return

    for query_name in query_names:
        if query_name not in lengths:
            print(f"Skipping '{query_name}': not found in {args.query_fasta}")
            continue
        query_rows = [r for r in rows if r["query_name"] == query_name]
        if not query_rows:
            print(f"Skipping '{query_name}': no hits in {args.csv}")
            continue

        hits = merge_regions_by_target(query_rows, args.gap_merge)
        n_regions = len(query_rows)
        n_targets = len({h["target_name"] for h in hits})
        if args.max_hits is not None and n_targets > args.max_hits:
            # Cap by distinct target, not by hit-span count: one heavily-fragmented
            # target (many small hits) would otherwise crowd out every other target.
            best_containment = {}
            for h in hits:
                best_containment[h["target_name"]] = max(
                    best_containment.get(h["target_name"], 0.0), h["containment"])
            top_targets = set(sorted(best_containment, key=best_containment.get,
                                      reverse=True)[: args.max_hits])
            hits = [h for h in hits if h["target_name"] in top_targets]

        base = os.path.join(args.output_dir, f"{safe_filename(query_name)}.hits")
        out_paths = [f"{base}.png", f"{base}.svg"]
        plot_gene(query_name, lengths[query_name], hits, out_paths, dpi=args.dpi)
        print(f"Wrote {base}.png / .svg "
              f"({len({h['target_name'] for h in hits})} targets, {len(hits)} hits, {n_regions} regions)")


if __name__ == "__main__":
    main()
