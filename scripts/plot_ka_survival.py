#!/usr/bin/env python3
"""Plot the distribution of the normalised region score x = lambda_pair * S that
`kmerseek index --ka-survival-out` wrote, with the fitted Karlin-Altschul line and the bins
it was read from.

    python scripts/plot_ka_survival.py survival_database.csv [survival_shuffled.csv ...] \\
        -o fig.png --subtitle "SCOPe40, hp_thomas_dill2 k = 12, C = 2, X = 8"

One column per CSV. Top row: regions in each bin of x, the points the line is fitted to.
Bottom row: the same regions summed (x at or above the bin). Grey band: the bins used.
Dotted line: 30 regions, the floor below which a bin is not fitted. The title of each
column gives the λ correction (1 means the closed-form per-pair lambda holds) with its
standard error over the fitted bins, and K.
"""

import argparse
import csv
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OBSERVED = "#1f77b4"
REFERENCE = "#7f7f7f"
FITTED = "#d62728"
WINDOW = "#cfcfcf"
FLOOR = "#111111"
MIN_BIN_COUNT = 30


def load(path):
    rows = list(csv.DictReader(open(path)))
    meta = rows[0]
    scores = [float(r["x"]) for r in rows]
    width = float(meta["bin_width"])
    survival = [int(r["n_regions_at_least"]) for r in rows]
    density = [survival[i] - (survival[i + 1] if i + 1 < len(survival) else 0) for i in range(len(survival))]
    window = [float(r["x"]) for r in rows if r["in_fit"] == "true"]
    ref_survival = [int(r["reference_n_regions_at_least"]) if r.get("reference_n_regions_at_least") else 0 for r in rows]
    ref_density = [ref_survival[i] - (ref_survival[i + 1] if i + 1 < len(ref_survival) else 0) for i in range(len(ref_survival))]
    # A refused fit (too few score bins above the peak) writes the curve with the fit
    # columns empty; the picture is then the histogram alone.
    if not meta["slope"]:
        return meta, scores, survival, density, window, None, None, None, ref_density, ref_survival, width
    lam, k = float(meta["slope"]), float(meta["k"])
    residues, kmers = float(meta["query_residues"]), float(meta["database_kmers"])
    per_bin = 1 - math.exp(-lam * width)
    ln_intercept = math.log(k * residues * kmers * per_bin)
    fit_density = [math.exp(ln_intercept - lam * s) for s in scores]
    fit_survival = [d / per_bin for d in fit_density]
    return meta, scores, survival, density, window, fit_density, fit_survival, ln_intercept, ref_density, ref_survival, width


def slope_standard_error(scores, density, window, lam, ln_intercept):
    points = [(s, d) for s, d in zip(scores, density) if s in window and d > 0]
    n = len(points)
    mean_s = sum(s for s, _ in points) / n
    sxx = sum((s - mean_s) ** 2 for s, _ in points)
    residuals = [math.log(d) - (ln_intercept - lam * s) for s, d in points]
    return math.sqrt(sum(r * r for r in residuals) / (n - 2) / sxx)


def draw_column(axes, path, first_column, args_label=None):
    meta, scores, survival, density, window, fit_density, fit_survival, ln_intercept, ref_density, ref_survival, width = load(path)
    fitted = fit_density is not None
    if fitted:
        lam, k = float(meta["slope"]), float(meta["k"])
        se = slope_standard_error(scores, density, window, lam, ln_intercept)
        xmax = max(window) + 12
    else:
        # No window to anchor on: show everything up to the last bin with a region.
        xmax = max(s for s, d in zip(scores, density) if d > 0)
    rows = [
        (density, fit_density, ref_density, "o", "", "regions in the bin (observed)"),
        (survival, fit_survival, ref_survival, "s", "none", "regions at or above the bin (observed, same regions summed)"),
    ]
    for ax, (observed, fit_line, reference, marker, face, label) in zip(axes, rows):
        if any(reference):
            keep = [i for i, o in enumerate(reference) if o > 0 and scores[i] <= xmax]
            ax.plot([scores[i] for i in keep], [reference[i] for i in keep], "^", ms=3.5, mfc="none",
                    color=REFERENCE, label="the same queries shuffled keeping dipeptides: the fit stops where the real curve rises above this")
        if fitted:
            ax.axvspan(min(window), max(window) + width, color=WINDOW, lw=0,
                       label=f"score bins the line is fitted on (≥ {MIN_BIN_COUNT} regions each)")
        keep = [i for i, o in enumerate(observed) if o > 0 and scores[i] <= xmax]
        style = dict(mfc=face) if face else {}
        ax.plot([scores[i] for i in keep], [observed[i] for i in keep], marker, ms=3.5,
                color=OBSERVED, label=label, **style)
        if fitted:
            keep = [i for i in range(len(scores)) if scores[i] <= xmax]
            ax.plot([scores[i] for i in keep], [fit_line[i] for i in keep], "-", color=FITTED, lw=1.8,
                    label="fitted line: its slope, sign flipped, is the λ correction; K from its height")
        ax.axhline(MIN_BIN_COUNT, color=FLOOR, ls=":", lw=1,
                   label=f"{MIN_BIN_COUNT} regions: bins below this are not fitted")
        ax.set_yscale("log")
        ax.set_ylim(0.5, max(observed) * 3)
        ax.spines[["top", "right"]].set_visible(False)
    label = args_label or f"{meta['null']} queries"
    verdict = (f"λ correction = {lam:.3f} ± {se:.3f}, K = {k:.4f}" if fitted
               else f"no fit: under 4 bins above the peak hold ≥ {MIN_BIN_COUNT} regions\n({meta['n_regions_at_least']} regions from {meta['n_queries']} quer{'y' if meta['n_queries'] == '1' else 'ies'})")
    axes[0].set_title(f"{label}\n{verdict}", fontsize=10, loc="left")
    axes[1].set_xlabel("x = λ_pair · S (nats)")
    if first_column:
        axes[0].set_ylabel(f"regions in each bin of x ({width:g} nat)\n(what the line is fitted to)")
        axes[1].set_ylabel("regions with x at or above the bin\n(the same fit, summed)")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", nargs="+")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--subtitle", default="")
    parser.add_argument("--labels", nargs="*", default=None, help="one title per CSV, in order")
    args = parser.parse_args()

    n = len(args.csv)
    fig, axes = plt.subplots(2, n, figsize=(4.5 * n, 7.6), dpi=150, sharex="col", squeeze=False)
    for j, path in enumerate(args.csv):
        label = args.labels[j] if args.labels and j < len(args.labels) else None
        draw_column([axes[0][j], axes[1][j]], path, first_column=(j == 0), args_label=label)
    handles = {}
    for ax in axes.flat:
        for h, l in zip(*ax.get_legend_handles_labels()):
            handles.setdefault(l, h)
    fig.legend(handles.values(), handles.keys(), loc="upper center", bbox_to_anchor=(0.5, 0.92),
               ncol=2, frameon=False, fontsize=9)
    title = "The λ correction is the slope and K the height of ln(regions at x); the fit stops where the real curve rises above the shuffled one"
    title += "\nS = matches − C × mismatches; λ_pair is the closed-form lambda for the pair's own compositions; ± is the slope's standard error"
    if args.subtitle:
        title += "\n" + args.subtitle
    fig.suptitle(title, fontsize=11, x=0.02, ha="left", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.83))
    fig.savefig(args.output, bbox_inches="tight")
    if args.output.endswith(".png"):
        fig.savefig(args.output[:-4] + ".svg", bbox_inches="tight")


if __name__ == "__main__":
    main()
