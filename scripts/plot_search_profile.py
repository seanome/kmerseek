"""Where one kmerseek search spends its time, step by step, on all of Swiss-Prot.

Reads docs/profiling/steps.tsv and timing.tsv, written by scripts/profile_search_steps.sbatch,
and draws one panel per index setting. Each step's seconds are the mean of the two runs with
the timers on.

    python scripts/plot_search_profile.py   # from the repository root
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import pubfig as pf  # noqa: E402

PROFILING = Path("docs/profiling")
OUT = Path("docs/images/search_step_profile_swissprot")
THREADS = 8

# Pipeline order, top to bottom, with the words used on the figure. Extension and
# Karlin-Altschul are left out: both are off by default, so neither ran (0 calls).
STEPS = {
    "LoadIndex": "Load index",
    "PrepareQuery": "Prepare query",
    "CandidateLookup": "Find candidate targets",
    "SignatureLoad": "Read target from disk",
    "Intersection": "Shared k-mers",
    "MatchedRegions": "Find matched regions",
    "RegionScoring": "Score regions",
    "ResultMetrics": "Other result columns",
    "SortResults": "Sort results",
    "WriteCsv": "Write CSV",
}
# Steps run on one thread; the rest run on 8 threads at once and are summed over them.
ONE_THREAD = {"LoadIndex", "WriteCsv"}
SETTINGS = [("protein20", 10), ("hp_lehninger2", 19)]


def mean_steps(steps, alphabet, k):
    """Mean seconds and calls per step over the profiled runs of one setting."""
    rows = steps[(steps.alphabet == alphabet) & (steps.k == k)]
    return rows.groupby("step")[["seconds", "calls"]].mean().reindex(list(STEPS))


def wall_seconds(timing, alphabet, k):
    """Mean wall time of the warm runs (run 1 reads the index from a cold disk cache)."""
    rows = timing[(timing.alphabet == alphabet) & (timing.k == k) & (timing.run > 1)]
    return rows.wall_s.mean(), rows.max_rss_kb.mean() * 1024 / 1e9


def draw_panel(ax, table, title):
    """Horizontal bars, one per step, with the seconds and the number of
    times the step ran (n) written at the end."""
    y = range(len(table))
    colours = [pf.OKABE_ITO["blue"] if s in ONE_THREAD else pf.OKABE_ITO["vermillion"]
               for s in table.index]
    ax.barh(list(y), table.seconds, color=colours, height=0.7)
    xmax = table.seconds.max()
    for yi, (step, row) in zip(y, table.iterrows()):
        seconds = f"{row.seconds:.2f}" if row.seconds < 1 else f"{row.seconds:,.0f}"
        ax.text(row.seconds + xmax * 0.02, yi, f"{seconds} s, n = {row.calls:,.0f}", va="center")
    ax.set_yticks(list(y), [STEPS[s] for s in table.index])
    ax.set_xlim(0, xmax * 1.55)
    ax.set_xlabel("Time (s)")
    ax.set_title(title, loc="left")


def main():
    pf.use_style()
    steps = pd.read_csv(PROFILING / "steps.tsv", sep="\t")
    timing = pd.read_csv(PROFILING / "timing.tsv", sep="\t")
    fig, axs = pf.figure(pf.TWO_COLUMN_MM, 70, ncols=2, sharey=True)
    for ax, (alphabet, k), letter in zip(axs, SETTINGS, "ab"):
        wall, rss_gb = wall_seconds(timing, alphabet, k)
        title = f"{alphabet}, k={k}: {wall:.0f} s wall, {rss_gb:.0f} GB peak memory"
        draw_panel(ax, mean_steps(steps, alphabet, k), title)
        pf.panel_label(ax, letter)
    axs[0].invert_yaxis()  # the y axis is shared, so turning it once turns both panels
    axs[0].barh([0], [0], color=pf.OKABE_ITO["blue"], label="one thread (wall time)")
    axs[0].barh([0], [0], color=pf.OKABE_ITO["vermillion"],
                label=f"{THREADS} threads, seconds summed over threads")
    pf.shared_legend(fig, ncol=2)
    pf.save(fig, str(OUT))


if __name__ == "__main__":
    main()
