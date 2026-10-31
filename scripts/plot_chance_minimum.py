"""Speed against pairs kept for each --max-shared-kmer-pvalue cutoff, hp_lehninger2 k=19 on
Swiss-Prot. Reads docs/profiling/chance_minimum/timing.tsv and kept_by_shared.tsv.

Panel a: wall time of the 100-query search at each cutoff. Panel b: the share of the
unfiltered pairs each cutoff keeps, for all pairs, for the pairs whose best region has the
largest score kmerseek can report, and for the 1% of pairs with the most shared k-mers in one
region.

    python scripts/plot_chance_minimum.py   # from the repository root
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import pubfig as pf  # noqa: E402

DATA = Path("docs/profiling/chance_minimum")
OUT = Path("docs/images/chance_minimum_swissprot_hp_lehninger2_k19")
CUTOFFS = {
    "none": "no cutoff",
    "1e-2": "0.01",
    "1e-4": "0.0001",
    "1e-6": "0.000001",
    "1e-8": "0.00000001",
}
GROUPS = {
    "pairs_kept": ("all pairs", pf.GREY, "o"),
    "max_score_kept": ("pairs with the largest region score", pf.OKABE_ITO["vermillion"], "s"),
    "top1pct_by_shared_kept": ("top 1% by shared k-mers in one region", pf.OKABE_ITO["blue"], "D"),
}


def load():
    timing = pd.read_csv(DATA / "timing.tsv", sep="\t", dtype={"cutoff": str})
    kept = pd.read_csv(DATA / "kept_by_shared.tsv", sep="\t", dtype={"cutoff": str})
    timing = timing[timing.alphabet == "hp_lehninger2"].set_index("cutoff").reindex(list(CUTOFFS))
    kept = kept[kept.alphabet == "hp_lehninger2"].set_index("cutoff").reindex(list(CUTOFFS))
    kept["pairs_kept"] = kept.pairs / kept.loc["none", "pairs"]
    return timing, kept


def draw_time(ax, timing):
    y = range(len(CUTOFFS))
    ax.barh(list(y), timing.wall_s, height=0.6, color=pf.GREY)
    for i, wall in enumerate(timing.wall_s):
        ax.text(wall + 5, i, f"{wall:.0f} s", va="center")
    ax.set_yticks(list(y), list(CUTOFFS.values()))
    ax.set_ylabel("Chance p-value cutoff")
    ax.set_xlim(0, 380)
    ax.set_xlabel("Wall time, 100 queries (s)")


def draw_kept(ax, kept):
    y = range(len(CUTOFFS))
    for column, (label, colour, marker) in GROUPS.items():
        ax.plot(kept[column] * 100, list(y), marker, color=colour, label=label)
    ax.set_xlim(0, 105)
    ax.set_xlabel("Unfiltered pairs kept (%)")
    ax.grid(axis="y")


def main():
    pf.use_style()
    timing, kept = load()
    fig, (a, b) = pf.figure(pf.TWO_COLUMN_MM, 55, ncols=2, sharey=True)
    draw_time(a, timing)
    draw_kept(b, kept)
    a.invert_yaxis()  # shared y axis: turning it once turns both panels
    pf.panel_label(a, "a")
    pf.panel_label(b, "b")
    pf.shared_legend(fig, ncol=3)
    pf.save(fig, str(OUT))


if __name__ == "__main__":
    main()
