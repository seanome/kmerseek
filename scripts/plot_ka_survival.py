#!/usr/bin/env python3
"""Plot the score distribution that `kmerseek index --ka-survival-out` wrote, with the
fitted Karlin-Altschul line and the score bins it was read from.

    python scripts/plot_ka_survival.py survival_database.csv [survival_shuffled.csv ...] \\
        -o fig.png --subtitle "SCOPe40, hp_thomas_dill2 k = 12, C = 2, X = 8"

One column per CSV. Top row: regions at each score S, the points the line is fitted to.
Bottom row: the same regions summed (score >= S). Grey band: the bins used. Dotted line:
30 regions, the floor below which a bin is not fitted. The title of each column gives
lambda with the standard error of the slope over the fitted bins.
"""

import argparse
import csv
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OBSERVED = "#1f77b4"
FITTED = "#d62728"
WINDOW = "#cfcfcf"
FLOOR = "#111111"
MIN_BIN_COUNT = 30


def load(path):
    rows = list(csv.DictReader(open(path)))
    meta = rows[0]
    scores = [int(r["score"]) for r in rows]
    survival = [int(r["n_regions_at_least"]) for r in rows]
    density = [survival[i] - (survival[i + 1] if i + 1 < len(survival) else 0) for i in range(len(survival))]
    window = [int(r["score"]) for r in rows if r["in_fit"] == "true"]
    lam, k = float(meta["lambda"]), float(meta["k"])
    residues, kmers = float(meta["query_residues"]), float(meta["database_kmers"])
    ln_intercept = math.log(k * residues * kmers * (1 - math.exp(-lam)))
    fit_density = [math.exp(ln_intercept - lam * s) for s in scores]
    fit_survival = [d / (1 - math.exp(-lam)) for d in fit_density]
    return meta, scores, survival, density, window, fit_density, fit_survival, ln_intercept


def slope_standard_error(scores, density, window, lam, ln_intercept):
    points = [(s, d) for s, d in zip(scores, density) if s in window and d > 0]
    n = len(points)
    mean_s = sum(s for s, _ in points) / n
    sxx = sum((s - mean_s) ** 2 for s, _ in points)
    residuals = [math.log(d) - (ln_intercept - lam * s) for s, d in points]
    return math.sqrt(sum(r * r for r in residuals) / (n - 2) / sxx)


def draw_column(axes, path, first_column):
    meta, scores, survival, density, window, fit_density, fit_survival, ln_intercept = load(path)
    lam, k = float(meta["lambda"]), float(meta["k"])
    se = slope_standard_error(scores, density, window, lam, ln_intercept)
    xmax = max(window) + 30
    rows = [
        (density, fit_density, "o", "", "regions with score exactly S (observed)"),
        (survival, fit_survival, "s", "none", "regions with score ≥ S (observed, same regions summed)"),
    ]
    for ax, (observed, fitted, marker, face, label) in zip(axes, rows):
        ax.axvspan(min(window) - 0.5, max(window) + 0.5, color=WINDOW, lw=0,
                   label=f"score bins the line is fitted on (≥ {MIN_BIN_COUNT} regions each)")
        keep = [i for i, o in enumerate(observed) if o > 0 and scores[i] <= xmax]
        style = dict(mfc=face) if face else {}
        ax.plot([scores[i] for i in keep], [observed[i] for i in keep], marker, ms=3.5,
                color=OBSERVED, label=label, **style)
        keep = [i for i in range(len(scores)) if scores[i] <= xmax]
        ax.plot([scores[i] for i in keep], [fitted[i] for i in keep], "-", color=FITTED, lw=1.8,
                label="fitted line, slope −λ")
        ax.axhline(MIN_BIN_COUNT, color=FLOOR, ls=":", lw=1,
                   label=f"{MIN_BIN_COUNT} regions: bins below this are not fitted")
        ax.set_yscale("log")
        ax.set_ylim(0.5, max(observed) * 3)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_title(f"{meta['null']} queries: λ = {lam:.3f} ± {se:.3f}, K = {k:.4f}", fontsize=10, loc="left")
    axes[1].set_xlabel(f"region score S = matches − {float(meta['mismatch_penalty']):g} × mismatches")
    if first_column:
        axes[0].set_ylabel("regions with score exactly S\n(what the line is fitted to)")
        axes[1].set_ylabel("regions with score ≥ S\n(the same fit, summed)")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", nargs="+")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--subtitle", default="")
    args = parser.parse_args()

    n = len(args.csv)
    fig, axes = plt.subplots(2, n, figsize=(4.5 * n, 7.6), dpi=150, sharex="col", squeeze=False)
    for j, path in enumerate(args.csv):
        draw_column([axes[0][j], axes[1][j]], path, first_column=(j == 0))
    handles = {}
    for ax in (axes[0][0], axes[1][0]):
        for h, l in zip(*ax.get_legend_handles_labels()):
            handles.setdefault(l, h)
    fig.legend(handles.values(), handles.keys(), loc="upper center", bbox_to_anchor=(0.5, 0.935),
               ncol=2, frameon=False, fontsize=9)
    title = "λ is the slope and K the intercept of ln(regions at score S); the fit stops where counts rise above the line"
    if args.subtitle:
        title += "\n" + args.subtitle + ", ± is the slope's standard error"
    fig.suptitle(title, fontsize=11, x=0.02, ha="left", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    fig.savefig(args.output, bbox_inches="tight")
    if args.output.endswith(".png"):
        fig.savefig(args.output[:-4] + ".svg", bbox_inches="tight")


if __name__ == "__main__":
    main()
