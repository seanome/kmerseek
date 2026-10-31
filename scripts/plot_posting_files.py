"""Swiss-Prot search before and after the inverted index moved to posting files read in place.

Reads docs/profiling/posting_files/timing.tsv and steps.tsv, written by
scripts/profile_posting_files.sbatch, and draws wall time, peak memory and the time to load
the index for each of the three runs of each version, one row per index setting. Run 1 read
the index from a cold disk cache; runs 2 and 3 found it already in memory.

    python scripts/plot_posting_files.py   # from the repository root
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import pubfig as pf  # noqa: E402

DATA = Path("docs/profiling/posting_files")
OUT = Path("docs/images/posting_files_swissprot")
SETTINGS = [("protein20", 10), ("hp_lehninger2", 19)]
# Grey for before and blue for after, as in docs/images/candidate_filter_*.png.
VERSIONS = {
    "before": ("before: HashMaps loaded at startup", pf.GREY, "o"),
    "after": ("after: posting files read in place", pf.OKABE_ITO["blue"], "s"),
}
RUN_LABELS = ["1\n(cold)", "2", "3"]


def per_run(timing, steps):
    """One row per (setting, version, run): wall seconds, peak GB and seconds to load the index."""
    keys = ["alphabet", "k", "version", "run"]
    load = steps[steps.step == "LoadIndex"][keys + ["seconds"]]
    table = timing.merge(load, on=keys)
    table["rss_gb"] = table.max_rss_kb * 1024 / 1e9
    return table.rename(columns={"seconds": "load_s"})


def number(value):
    """One decimal below 10, none above, so values near zero still read."""
    return f"{value:.1f}" if value < 10 else f"{value:.0f}"


def draw_panel(ax, table, column, ylabel):
    """Each version's value per run as a point, its number written above it."""
    for offset, (version, (label, colour, marker)) in zip((-0.12, 0.12), VERSIONS.items()):
        rows = table[table.version == version].sort_values("run")
        x = rows.run - 1 + offset
        ax.plot(x, rows[column], marker, color=colour, label=label, linestyle="none")
        for xi, value in zip(x, rows[column]):
            ax.annotate(number(value), (xi, value), xytext=(0, 4),
                        textcoords="offset points", ha="center", va="bottom")
    ax.set_xticks(range(len(RUN_LABELS)), RUN_LABELS)
    ax.set_xlim(-0.5, len(RUN_LABELS) - 0.3)
    ax.set_ylim(0, table[column].max() * 1.2)
    ax.set_xlabel("Run")
    ax.set_ylabel(ylabel)


def main():
    pf.use_style()
    timing = pd.read_csv(DATA / "timing.tsv", sep="\t")
    steps = pd.read_csv(DATA / "steps.tsv", sep="\t")
    table = per_run(timing, steps)
    fig, axs = pf.figure(pf.TWO_COLUMN_MM, 110, nrows=2, ncols=3)
    panels = [("wall_s", "Search wall time (s)"),
              ("rss_gb", "Peak memory (GB)"),
              ("load_s", "Time to load the index (s)")]
    letters = iter("abcdef")
    for row, (alphabet, k) in zip(axs, SETTINGS):
        rows = table[(table.alphabet == alphabet) & (table.k == k)]
        for ax, (column, ylabel) in zip(row, panels):
            draw_panel(ax, rows, column, ylabel)
            pf.panel_label(ax, next(letters))
        row[0].set_title(f"{alphabet}, k={k}", loc="left", pad=12)
    pf.shared_legend(fig, ncol=2)
    pf.save(fig, str(OUT))


if __name__ == "__main__":
    main()
