"""What filtering candidates before reading targets saves, hp_lehninger2 k=19 on Swiss-Prot.

Panel a: thread-seconds per step before (docs/profiling/steps.tsv, runs 2 and 4) and after
(docs/profiling/candidate_filter/steps.tsv, uncapped run 2). Panel b: wall time for each
--max-kmer-targets cap, and the share of uncapped query-target pairs it keeps.

    python scripts/plot_candidate_filter.py   # from the repository root
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import pubfig as pf  # noqa: E402

BEFORE = Path("docs/profiling")
AFTER = BEFORE / "candidate_filter"
OUT = Path("docs/images/candidate_filter_swissprot_hp_lehninger2_k19")
STEPS = {
    "SignatureLoad": "Read target from disk",
    "Intersection": "Shared k-mers",
    "MatchedRegions": "Find matched regions",
}
CAP_LABELS = {"none": "no cap", "10000": "10,000", "1000": "1,000", "100": "100"}


def step_seconds():
    """Thread-seconds per step, before (mean of the profiled runs) and after (warm run)."""
    before = pd.read_csv(BEFORE / "steps.tsv", sep="\t")
    before = before[before.alphabet == "hp_lehninger2"].groupby("step").seconds.mean()
    after = pd.read_csv(AFTER / "steps.tsv", sep="\t")
    after = after[(after.alphabet == "hp_lehninger2") & (after.cap == "none") & (after.run == 2)]
    after = after.set_index("step").seconds
    return before.reindex(list(STEPS)), after.reindex(list(STEPS))


def draw_steps(ax):
    before, after = step_seconds()
    y = range(len(STEPS))
    ax.barh([i - 0.18 for i in y], before, height=0.36, color=pf.GREY, label="before")
    ax.barh([i + 0.18 for i in y], after, height=0.36, color=pf.OKABE_ITO["blue"],
            label="count shared k-mers before reading targets")
    for i, (b, a) in enumerate(zip(before, after)):
        ax.text(b + 10, i - 0.18, f"{b:,.0f} s", va="center")
        ax.text(a + 10, i + 0.18, f"{a:,.0f} s", va="center")
    ax.set_yticks(list(y), list(STEPS.values()))
    ax.invert_yaxis()
    ax.set_xlim(0, 760)
    ax.set_xlabel("Time summed over 8 threads (s)")


def draw_caps(ax):
    """One bar per cap: wall time, with the share of uncapped pairs it keeps written beside."""
    timing = pd.read_csv(AFTER / "timing.tsv", sep="\t", dtype={"cap": str})
    kept = pd.read_csv(AFTER / "kept.tsv", sep="\t", dtype={"cap": str}).set_index("cap")
    # Run 1 without a cap read the index from a cold disk cache; leave it out.
    hp = timing[(timing.alphabet == "hp_lehninger2") & ~((timing.cap == "none") & (timing.run == 1))]
    wall = hp.groupby("cap").wall_s.mean().reindex(list(CAP_LABELS))
    y = range(len(CAP_LABELS))
    ax.barh(list(y), wall, height=0.6, color=pf.OKABE_ITO["blue"])
    for i, cap in enumerate(CAP_LABELS):
        share = kept.loc[cap, "pairs_kept_fraction"] * 100
        ax.text(wall[cap] + 4, i, f"{wall[cap]:.0f} s, {share:.1f}% of pairs kept", va="center")
    ax.set_yticks(list(y), [f"cap {v}" if k != "none" else v for k, v in CAP_LABELS.items()])
    ax.invert_yaxis()
    ax.set_xlim(0, 330)
    ax.set_xlabel("Wall time (s)")


def main():
    pf.use_style()
    fig, (a, b) = pf.figure(pf.TWO_COLUMN_MM, 60, ncols=2, width_ratios=[1.3, 1])
    draw_steps(a)
    draw_caps(b)
    pf.panel_label(a, "a")
    pf.panel_label(b, "b")
    pf.shared_legend(fig, ncol=2)
    pf.save(fig, str(OUT))


if __name__ == "__main__":
    main()
