"""Plot bench/ledger.tsv: search time, search memory and index size per tool.

One bar per tool at the mean of the ledger rows for the given commits, with each row
drawn as a dot on it, so run-to-run spread shows.

    python bench/plot_ledger.py bench/ledger.tsv docs/images/bench_speed_baseline bc3348a f367d3c
"""

import sys
from pathlib import Path

import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import pubfig as pf  # noqa: E402

TOOLS = {"kmerseek": "kmerseek", "diamond": "DIAMOND", "mmseqs": "MMseqs2"}
COLOURS = {
    "kmerseek": pf.OKABE_ITO["blue"],
    "diamond": pf.OKABE_ITO["vermillion"],
    "mmseqs": pf.OKABE_ITO["bluish_green"],
}
SEARCH_STEP = {"kmerseek": "search", "diamond": "blastp", "mmseqs": "search"}
INDEX_STEP = {"kmerseek": "index", "diamond": "makedb", "mmseqs": "createindex"}


def panel_values(ledger):
    """Per tool: the search rows' wall time and peak memory, and index size in GB."""
    rows = []
    for tool in TOOLS:
        mine = ledger[ledger.tool == tool]
        search = mine[mine.step == SEARCH_STEP[tool]]
        index = mine[mine.step == INDEX_STEP[tool]]
        rows.append(
            {
                "tool": tool,
                "Search wall time (s)": list(search.wall_s),
                "Search peak memory (GB)": list(search.peak_rss_gb),
                "Index on disk (GB)": list(index.index_bytes.astype(float) / 1024**3),
            }
        )
    return pd.DataFrame(rows)


def draw(ax, values, column):
    for i, row in values.iterrows():
        ys = row[column]
        ax.bar(i, sum(ys) / len(ys), color=COLOURS[row.tool], width=0.6)
        ax.plot([i] * len(ys), ys, "o", color="black", markersize=2)
    ax.set_xticks(range(len(values)), [TOOLS[t] for t in values.tool])
    ax.set_ylabel(column)


def main(ledger_path, out_stem, *commits):
    ledger = pd.read_csv(ledger_path, sep="\t", dtype={"git_commit": str})
    if commits:
        ledger = ledger[ledger.git_commit.isin(commits)]
    values = panel_values(ledger)
    pf.use_style()
    fig, axs = pf.figure(pf.TWO_COLUMN_MM, 55, ncols=3)
    columns = ["Search wall time (s)", "Search peak memory (GB)", "Index on disk (GB)"]
    key = [
        Patch(color=pf.GREY, label="bar: mean of the runs"),
        Line2D(
            [],
            [],
            marker="o",
            color="black",
            linestyle="",
            markersize=2,
            label="dot: one run (one ledger row)",
        ),
    ]
    fig.legend(handles=key, loc="outside upper center", ncol=2)
    for ax, column, letter in zip(axs, columns, "abc"):
        draw(ax, values, column)
        pf.panel_label(ax, letter)
    pf.save(fig, out_stem)


if __name__ == "__main__":
    main(*sys.argv[1:])
