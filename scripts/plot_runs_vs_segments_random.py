#!/usr/bin/env python3
"""Why one exact region gets two different E-values: on a random diagonal, chance makes
far more +1/-2 regions scoring at least S than exact runs of at least S matches.

    python scripts/plot_runs_vs_segments_random.py -o docs/images/runs_vs_segments_random_pr_same_half

Each position of one long diagonal is in the same class for query and target with
probability 0.5, independently (hp_lehninger2 pairs have Pr(same) 0.498 to 0.502).
An exact run is counted once, where it starts: what region_run_evalue counts. A +1/-2
region is counted once per climb of the running score, which restarts at 0: what
region_ka_evalue counts with --extend-mismatch-penalty 2.
"""

import argparse
import random

import matplotlib

matplotlib.use("Agg")
import pubfig as pf

POSITIONS = 20_000_000
PR_SAME = 0.5
PENALTY = 2
SCORES = range(4, 21, 2)


def count_regions(seed=1):
    """Best score of every exact run and of every +1/-2 climb along one random diagonal."""
    rng = random.Random(seed)
    runs, climbs = [], []
    run = score = best = 0
    for _ in range(POSITIONS):
        if rng.random() < PR_SAME:
            run += 1
            score += 1
            best = max(best, score)
            continue
        if run:
            runs.append(run)
        run = 0
        score -= PENALTY
        if score <= 0:
            if best:
                climbs.append(best)
            score = best = 0
    return runs, climbs


def per_million(values, s):
    return 1e6 * sum(v >= s for v in values) / POSITIONS


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("-o", "--output", required=True, help="path without extension")
    args = parser.parse_args()
    runs, climbs = count_regions()
    pf.use_style()
    fig, ax = pf.figure(pf.ONE_COLUMN_MM, 70)
    ax.set_yscale("log")
    formula = [1e6 * (1 - PR_SAME) * PR_SAME**s for s in SCORES]
    ax.plot(
        SCORES, formula, color=pf.GREY, label="(1 - Pr(same)) x Pr(same)^S, the run E-value formula"
    )
    ax.plot(
        SCORES,
        [per_million(runs, s) for s in SCORES],
        "o",
        color=pf.OKABE_ITO["blue"],
        label="exact runs of at least S matches (region_run_evalue counts these)",
    )
    ax.plot(
        SCORES,
        [per_million(climbs, s) for s in SCORES],
        "^",
        color=pf.OKABE_ITO["vermillion"],
        label="regions scoring at least S at +1 / -2 (region_ka_evalue counts these)",
    )
    ax.set_xticks(list(SCORES))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel("S: run length, or region score at +1 per match, -2 per mismatch")
    ax.set_ylabel("Found by chance per million positions\n(random sequences, Pr(same) = 0.5)")
    pf.shared_legend(fig, ncol=1)
    pf.save(fig, args.output, formats=("png", "svg"))


if __name__ == "__main__":
    main()
