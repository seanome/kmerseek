#!/usr/bin/env python3
"""Plot the score distribution that `kmerseek index --ka-survival-out` wrote, with the
fitted Karlin-Altschul line and the score range it was read from.

    python scripts/plot_ka_survival.py survival_database.csv [survival_shuffled.csv ...] -o fig.png

One panel per CSV. Dots: regions with score >= S. Line: the fit, K L N e^(-lambda S).
Shaded band: the score bins the line was fitted on. Dashed vertical line: the first score
bin above the fit whose count sat above the line, where related sequences start to show.
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
EXCESS = "#111111"


def read_curve(path):
    rows = list(csv.DictReader(open(path)))
    meta = rows[0]
    scores = [int(r["score"]) for r in rows]
    observed = [int(r["n_regions_at_least"]) for r in rows]
    fitted = [float(r["fitted_n_regions_at_least"]) for r in rows]
    window = [int(r["score"]) for r in rows if r["in_fit"] == "true"]
    return meta, scores, observed, fitted, window


def excess_start(scores, observed, fitted, window):
    """First score above the fit window where the observed count sits above the line by
    more than two Poisson standard deviations, or None."""
    for s, o, f in zip(scores, observed, fitted):
        if s <= max(window) or o < 1:
            continue
        if math.log(o) - math.log(f) > 2.0 / math.sqrt(o) + 0.05:
            return s
    return None


def draw_panel(ax, path, show_ylabel):
    meta, scores, observed, fitted, window = read_curve(path)
    lam = float(meta["lambda"])
    k = float(meta["k"])
    penalty = float(meta["mismatch_penalty"])
    keep = [i for i, o in enumerate(observed) if o >= 1]
    ax.axvspan(min(window) - 0.5, max(window) + 0.5, color=WINDOW, lw=0,
               label="score bins the line was fitted on")
    ax.plot([scores[i] for i in keep], [observed[i] for i in keep], "o", ms=3.5,
            color=OBSERVED, label="regions with score ≥ S (observed)")
    ax.plot(scores, fitted, "-", color=FITTED, lw=1.8,
            label="fitted line K·L·N·e^(−λS)")
    start = excess_start(scores, observed, fitted, window)
    if start is not None:
        ax.axvline(start, color=EXCESS, ls="--", lw=1.2,
                   label="first score where counts rise above the line (related sequences)")
    ax.set_yscale("log")
    ax.set_ylim(0.5, max(observed) * 3)
    ax.set_xlim(min(scores) - 1, max(min(scores) + 60, max(window) + 15))
    ax.set_title(f"{meta['null']} queries: λ = {lam:.3f}, K = {k:.4f}", fontsize=10, loc="left")
    ax.set_xlabel(f"region score S = matches − {penalty:g} × mismatches")
    if show_ylabel:
        ax.set_ylabel(f"regions with score ≥ S\n({meta['n_queries']} calibration queries)")
    ax.spines[["top", "right"]].set_visible(False)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", nargs="+")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--title", default="λ and K are read off the straight part of the curve; related sequences lift the right-hand end")
    args = parser.parse_args()

    n = len(args.csv)
    fig, axes = plt.subplots(1, n, figsize=(4.6 * n, 4.4), dpi=150, sharey=True, squeeze=False)
    for i, (ax, path) in enumerate(zip(axes[0], args.csv)):
        draw_panel(ax, path, show_ylabel=(i == 0))
    handles, labels = axes[0][0].get_legend_handles_labels()
    seen = {}
    for h, l in zip(handles, labels):
        seen.setdefault(l, h)
    fig.legend(seen.values(), seen.keys(), loc="upper center", bbox_to_anchor=(0.5, 0.93),
               ncol=2, frameon=False, fontsize=9)
    fig.suptitle(args.title, fontsize=11, x=0.02, ha="left", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.84))
    fig.savefig(args.output, bbox_inches="tight")
    if args.output.endswith(".png"):
        fig.savefig(args.output[:-4] + ".svg", bbox_inches="tight")


if __name__ == "__main__":
    main()
