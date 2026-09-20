#!/usr/bin/env python3
"""Draw every matched region of a query/target pair as a dot plot, next to a histogram of
the regions' p-values against the default region cutoff (p = 0.05).

Input is the CSV from `kmerseek search -o hits.csv` (one row per matched region). A second
CSV from an older build overlays its p-values as an outline, to show what a change to region
detection did to one pair. The figure in docs/images was made with

    kmerseek index -i tests/testdata/fasta/bcl2_first25_*.fasta.gz -o idx -k 9 -a hp_lehninger2
    kmerseek search -q tests/testdata/fasta/bcl2_first25_*.fasta.gz -t idx -k 9 -a hp_lehninger2 -o after.csv
    # same search with the build from the parent commit -> before.csv
    python scripts/plot_matched_regions.py --after after.csv --before before.csv \\
        --pair BCL2_HUMAN RTN3_HUMAN --pair B2LA1_HUMAN ASPP2_HUMAN \\
        --query-fasta tests/testdata/fasta/bcl2_first25_*.fasta.gz \\
        --output docs/images/regions_before_after_chaining
"""

import argparse
import gzip

import matplotlib
import numpy as np
import polars as pl

matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

REGION = "#8a8a8a"  # a matched region in the current build
BEFORE = "#e08214"  # regions from the older build, histogram outline only
BEST = "#2166ac"  # the best-scoring region (lowest p) in the current build
CUTOFF = 0.05
HIST_BINS = np.logspace(np.log10(0.03), 0, 25)


def read_fasta_lengths(path):
    """Return {short name: length in aa}, keyed by the UniProt entry name (BCL2_HUMAN)."""
    opener = gzip.open if path.endswith(".gz") else open
    lengths, name = {}, None
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0].split("|")[-1]
                lengths[name] = 0
            else:
                lengths[name] += len(line.strip())
    return lengths


def pair_rows(df, query, target):
    return df.filter(
        pl.col("query_name").str.contains(query) & pl.col("target_name").str.contains(target)
    )


def draw_dot_plot(ax, rows, best, query, qlen, target, tlen, n_before):
    for r in rows.iter_rows(named=True):
        ax.plot(
            [r["region_start"], r["region_end"]],
            [r["target_start"], r["target_end"]],
            color=REGION,
            lw=2,
            solid_capstyle="butt",
        )
    ax.plot(
        [best["region_start"], best["region_end"]],
        [best["target_start"], best["target_end"]],
        color=BEST,
        lw=4,
        solid_capstyle="butt",
    )
    ax.set_xlim(0, qlen)
    ax.set_ylim(0, tlen)
    ax.set_xlabel(f"position in query {query.split('_')[0]} (aa)")
    ax.set_ylabel(f"position in target {target.split('_')[0]} (aa)")
    ax.set_title(
        f"{query.split('_')[0]} ({qlen} aa) vs {target.split('_')[0]} ({tlen} aa)\n"
        f"{rows['n_intersecting_hashes'][0]} shared k-mers in {rows.height} regions"
        f" (were {n_before} in the older build)",
        fontsize=10,
        loc="left",
    )
    ax.annotate(
        f"best region: {best['region_length']} aa, p = {best['region_tail_probability']:.3f}",
        xy=(best["region_end"], best["target_end"]),
        xytext=(0.97, 0.06),
        textcoords="axes fraction",
        ha="right",
        va="bottom",
        color=BEST,
        fontsize=9,
        bbox=dict(boxstyle="square,pad=0.2", facecolor="white", edgecolor="none", alpha=0.9),
        arrowprops=dict(arrowstyle="-", color=BEST, lw=0.8),
    )


def draw_pvalue_histogram(ax, rows, before, best):
    ax.hist(rows["region_tail_probability"], bins=HIST_BINS, color=REGION)
    ax.hist(before["region_tail_probability"], bins=HIST_BINS, histtype="step", edgecolor=BEFORE, lw=1.5)
    ax.axvline(CUTOFF, color="black", ls="--", lw=1)
    best_p = best["region_tail_probability"]
    ax.axvline(best_p, color=BEST, lw=3)
    ax.set_xscale("log")
    ticks = [0.03, 0.05, 0.1, 0.2, 0.5, 1.0]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_xlabel(
        "p-value of the region\n"
        "(Poisson tail: chance of at least this many shared k-mers in a window this long)",
        fontsize=8.5,
    )
    ax.set_ylabel("number of regions")
    verdict = "passes" if best_p < CUTOFF else "does not pass"
    ax.set_title(f"best region p = {best_p:.4f}, cutoff {CUTOFF}: {verdict}", fontsize=10, loc="left")


def legend_handles():
    return [
        Line2D([], [], color=REGION, lw=2, label="a matched region in this build: shared k-mers on one diagonal"),
        Patch(facecolor="none", edgecolor=BEFORE, lw=1.5, label="regions in the older build"),
        Line2D([], [], color=BEST, lw=3, label="best-scoring region (lowest p) in this build"),
        Line2D([], [], color="black", ls="--", lw=1, label=f"default region cutoff, p = {CUTOFF}"),
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--after", required=True, help="hits CSV from the current build")
    parser.add_argument("--before", required=True, help="hits CSV from the older build")
    parser.add_argument("--pair", nargs=2, action="append", metavar=("QUERY", "TARGET"), required=True)
    parser.add_argument("--query-fasta", required=True, help="FASTA with both proteins, for their lengths")
    parser.add_argument("--title", default="")
    parser.add_argument("--output", required=True, help="path without extension; writes .png and .svg")
    args = parser.parse_args()

    after, before = pl.read_csv(args.after), pl.read_csv(args.before)
    lengths = read_fasta_lengths(args.query_fasta)
    n = len(args.pair)
    fig, axes = plt.subplots(
        n, 2, figsize=(11, 4.75 * n), squeeze=False,
        gridspec_kw={"width_ratios": [1, 1.1], "hspace": 0.5, "wspace": 0.3},
    )
    if args.title:
        fig.suptitle(args.title, fontsize=11, y=0.985)
    fig.legend(handles=legend_handles(), loc="upper center", bbox_to_anchor=(0.5, 0.905), ncol=2, fontsize=9, frameon=False)

    for row, (query, target) in enumerate(args.pair):
        rows, old = pair_rows(after, query, target), pair_rows(before, query, target)
        best = rows.sort("region_tail_probability").row(0, named=True)
        draw_dot_plot(axes[row, 0], rows, best, query, lengths[query], target, lengths[target], old.height)
        draw_pvalue_histogram(axes[row, 1], rows, old, best)
        for ax in axes[row]:
            ax.spines[["top", "right"]].set_visible(False)

    fig.subplots_adjust(top=0.82, bottom=0.07, left=0.08, right=0.98)
    fig.savefig(f"{args.output}.png", dpi=150)
    fig.savefig(f"{args.output}.svg")


if __name__ == "__main__":
    main()
