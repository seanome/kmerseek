"""Draw the scores and fits of the tests in src/rust/evalue.rs, one figure per test.

`fit_scores` reads lambda and K off the regions of a search: it counts the regions in
each score bin, and fits a line to ln(count) against score over the bins above the most
populated one that hold at least 30 regions. It grows the fit one bin at a time and stops
at the first bin that sits too far above the line (where related pairs begin).
`fit_scores_with_reference` does the same on the ratio of each bin's count to the count of
a reference: the same queries, shuffled.

The tests write their counts and the fit they got as JSON when KMERSEEK_FIGURE_DATA is
set, so the figures show what the Rust code returned, not a copy of it:

    mkdir -p /tmp/fit_tests
    KMERSEEK_FIGURE_DATA=/tmp/fit_tests cargo test --no-default-features --lib evalue::
    python scripts/plot_fit_scores_tests.py /tmp/fit_tests --output docs/images/fit_scores_tests
"""

import argparse
import json
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, NullLocator

COUNT_COLOR = "#444444"  # this test's regions per bin
FIT_COLOR = "#2a9d8f"  # the bins the line was fitted on, and the line
BEND_COLOR = "#e76f51"  # the bin where the fit stopped
REFERENCE_COLOR = "#7b5ea7"  # the reference (shuffled) regions per bin
FLOOR_COLOR = "#9e9e9e"  # the fewest regions a bin needs
SURVIVAL_COLOR = "#8ab6d6"  # regions at or above a bin


def usable_bins(bins, reference_bins, min_count):
    """Bins above the most populated one with at least `min_count` regions, as in
    `tail_bins` plus the count filter of `fit_scores` (and, given a reference, of
    `fit_scores_with_reference`). Returns (score, count, reference count or None)."""
    if not bins:
        return []
    peak = max(range(len(bins)), key=lambda i: (bins[i][1], -i))
    reference = dict(map(tuple, reference_bins))
    out = []
    for score, count in bins[peak + 1 :]:
        ref = reference.get(score) if reference else None
        if count < min_count or (reference and (ref is None or ref < min_count)):
            continue
        out.append((score, count, ref))
    return out


def least_squares(points):
    """Slope and intercept of y on x."""
    n = len(points)
    mx = sum(x for x, _ in points) / n
    my = sum(y for _, y in points) / n
    sxx = sum((x - mx) ** 2 for x, _ in points)
    sxy = sum((x - mx) * (y - my) for x, y in points)
    slope = sxy / sxx if sxx > 0 else 0.0
    return slope, my - slope * mx


def ratio_baseline(usable, base):
    """The line `fit_scores_with_reference` measures the ratio against: least squares of
    ln(count / reference count) on score over the lowest `base` usable bins."""
    points = [(s, math.log(c) - math.log(r)) for s, c, r in usable[:base]]
    return least_squares(points)


def fit_label(fit, made_with):
    text = f"fitted lambda {fit['lambda']:.3f} per bin"
    if made_with is not None:
        text += f" (counts made with {made_with:g})"
    return text


def draw_counts(ax, record, case):
    """Regions per score bin on a log axis, with the bins the fit used and the line."""
    fit = case["fit"]
    window = range(fit["score_lo"], fit["score_hi"] + 1) if fit else range(0)
    ax.axhline(record["min_bin_count"], color=FLOOR_COLOR, ls=":", lw=1.2)
    if case["reference_bins"]:
        xs, ys = zip(*[(s, c) for s, c in case["reference_bins"] if c > 0])
        ax.plot(xs, ys, "s", ms=4, mfc="none", color=REFERENCE_COLOR)
    for score, count in case["bins"]:
        if count == 0:
            continue
        inside = score in window
        ax.plot(
            score,
            count,
            "o",
            ms=6,
            color=FIT_COLOR if inside else COUNT_COLOR,
            mfc=FIT_COLOR if inside else "none",
        )
    if fit:
        xs = [fit["score_lo"], fit["score_hi"]]
        ax.plot(
            xs,
            [math.exp(fit["ln_intercept"] - fit["lambda"] * x) for x in xs],
            color=FIT_COLOR,
            lw=1.5,
        )
        if fit["bend_score"] is not None:
            bend = dict(map(tuple, case["bins"]))[fit["bend_score"]]
            ax.plot(fit["bend_score"], bend, "^", ms=11, color=BEND_COLOR, mfc="none", mew=2)
    ax.set_yscale("log")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("score bin")
    ax.set_ylabel("regions in the bin (log scale)")


def draw_ratio(ax, record, case):
    """ln(count / reference count) per usable bin, the baseline, and how far above it a
    bin may sit and still join the fit."""
    usable = usable_bins(case["bins"], case["reference_bins"], record["min_bin_count"])
    slope, intercept = ratio_baseline(usable, record["reference_base"])
    fit = case["fit"]
    base_scores = [s for s, _, _ in usable[: record["reference_base"]]]
    xs = [usable[0][0], usable[-1][0]]
    ax.plot(xs, [intercept + slope * x for x in xs], color=REFERENCE_COLOR, lw=1.2, ls="--")
    for score, count, ref in usable:
        y = math.log(count) - math.log(ref)
        allowed = intercept + slope * score
        allowed += record["bend_sigmas"] * math.sqrt(1 / count + 1 / ref) + record["bend_slack"]
        ax.plot([score - 0.35, score + 0.35], [allowed] * 2, color=FLOOR_COLOR, lw=1.5)
        inside = fit and fit["score_lo"] <= score <= fit["score_hi"]
        ax.plot(
            score,
            y,
            "o",
            ms=6,
            color=FIT_COLOR if inside else COUNT_COLOR,
            mfc=FIT_COLOR if inside else "none",
        )
        if fit and fit["bend_score"] == score:
            ax.plot(score, y, "^", ms=11, color=BEND_COLOR, mfc="none", mew=2)
    ax.axvspan(base_scores[0] - 0.5, base_scores[-1] + 0.5, color=REFERENCE_COLOR, alpha=0.08)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("score bin")
    ax.set_ylabel("ln(regions / reference regions)")


def plain_count_ticks(ax):
    """Label a log count axis spanning 20 to 200 regions with plain numbers."""
    ax.set_ylim(20, 200)
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_yticks([20, 30, 40, 60, 100, 200], labels=["20", "30", "40", "60", "100", "200"])


def count_legend(reference=False):
    handles = [
        Line2D(
            [],
            [],
            marker="o",
            ls="",
            color=COUNT_COLOR,
            mfc="none",
            label="regions in a bin the line does not use",
        ),
        Line2D(
            [],
            [],
            marker="o",
            ls="",
            color=FIT_COLOR,
            label="regions in a bin the line was fitted on",
        ),
        Line2D([], [], color=FIT_COLOR, lw=1.5, label="fitted line, slope -lambda"),
        Line2D(
            [],
            [],
            marker="^",
            ls="",
            ms=10,
            color=BEND_COLOR,
            mfc="none",
            mew=2,
            label="first bin too far above the line: fit stops",
        ),
        Line2D([], [], color=FLOOR_COLOR, ls=":", label="fewest regions a bin needs (30)"),
    ]
    if reference:
        handles.append(
            Line2D(
                [],
                [],
                marker="s",
                ls="",
                color=REFERENCE_COLOR,
                mfc="none",
                label="reference (shuffled) regions in the bin",
            )
        )
    return handles


def ratio_legend():
    return [
        Line2D(
            [],
            [],
            color=REFERENCE_COLOR,
            ls="--",
            label="baseline: line through the 8 lowest usable bins (shaded)",
        ),
        Line2D(
            [],
            [],
            color=FLOOR_COLOR,
            lw=1.5,
            label="highest ratio a bin may have and join:\nbaseline + 2 sqrt(1 / regions + 1 / reference regions) + 0.05",
        ),
    ]


def figure_with_legend(nrows, ncols, handles, title, width=5.2, height=3.8):
    """A grid of panels under a title and a two-column legend, the legend above the marks."""
    title_inches = 0.1 + 0.25 * (title.count("\n") + 1)
    legend_inches = title_inches + 0.1 + 0.17 * math.ceil(len(handles) / 2)
    fig_height = height * nrows + legend_inches
    fig, axes = plt.subplots(nrows, ncols, figsize=(width * ncols, fig_height), squeeze=False)
    fig.suptitle(title, y=1 - 0.08 / fig_height, va="top", fontsize=12)
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1 - title_inches / fig_height),
        ncol=2,
        fontsize=8.5,
        frameon=False,
    )
    # Read by main(), which lays the panels out once their titles are set.
    fig.panels_top = 1 - legend_inches / fig_height
    return fig, axes


def plot_counts_test(record):
    """`test_survival_counts` and `test_bin_counts`: the two ways of counting regions."""
    handles = [
        Line2D(
            [],
            [],
            marker="s",
            ls="",
            ms=9,
            color=COUNT_COLOR,
            label="bin_counts: regions with score in the bin",
        ),
        Line2D(
            [],
            [],
            marker="s",
            ls="",
            ms=9,
            color=SURVIVAL_COLOR,
            label="survival_counts: regions with score >= the bin",
        ),
    ]
    fig, axes = figure_with_legend(
        1,
        2,
        handles,
        "Scores fall into the bin below them; infinite and undefined scores are left out",
        height=3.2,
    )
    for ax, case in zip(axes[0], record["cases"]):
        survival, running = [], 0
        for score, count in reversed(case["bins"]):
            running += count
            survival.append((score, running))
        survival = dict(survival)
        for score, count in case["bins"]:
            for dx, n, color in (
                (-0.18, count, COUNT_COLOR),
                (0.18, survival[score], SURVIVAL_COLOR),
            ):
                ax.bar(score + dx, n, width=0.34, color=color)
                ax.text(score + dx, n + 0.08, str(n), ha="center", va="bottom", fontsize=9)
        ax.set_title(f"scores: {case['label']}", fontsize=10)
        ax.set_xticks([s for s, _ in case["bins"]])
        ax.set_xlabel("score bin")
        ax.set_ylabel("regions")
        ax.set_ylim(0, max(survival.values()) + 1)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    return fig


def plot_single_curve_test(record):
    """`test_fit_scores_recovers_slope_and_stops_at_homolog_excess`."""
    clean, bent = record["cases"]
    title = (
        f"The line fits bins {bent['fit']['score_lo']}-{bent['fit']['score_hi']} with or without "
        f"400 extra regions at score {bent['fit']['bend_score']}, where the fit stops"
    )
    fig, axes = figure_with_legend(1, 2, count_legend(), title)
    for ax, case in zip(axes[0], (clean, bent)):
        draw_counts(ax, record, case)
        ax.set_title(
            f"{case['label']}\n{fit_label(case['fit'], case['made_with_lambda'])}", fontsize=10
        )
    return fig


def plot_reference_test(record):
    """`test_fit_scores_with_reference_stops_where_the_ratio_rises`."""
    ramped, plain = record["cases"]
    title = (
        f"Related pairs rising slowly from score 24 stop the ratio test at {ramped['fit']['bend_score']}; "
        "without them the fit runs until bins hold fewer than 30 regions"
    )
    fig, axes = figure_with_legend(
        2, 2, count_legend(reference=True) + ratio_legend(), title, width=6.0, height=4.0
    )
    for row, case in zip(axes, (ramped, plain)):
        draw_counts(row[0], record, case)
        fit = case["fit"]
        row[0].set_title(
            f"{case['label']}\n{fit_label(fit, case['made_with_lambda'])}; reference {fit['reference_lambda']:.3f}",
            fontsize=10,
        )
        draw_ratio(row[1], record, case)
        row[1].set_title(f"{case['label']}: ratio to the reference", fontsize=10)
    return fig


def plot_too_few_bins_test(record):
    """`test_fit_scores_needs_enough_bins`."""
    case = record["cases"][0]
    fig, axes = figure_with_legend(
        1,
        1,
        count_legend()[:1] + count_legend()[4:],
        "3 bins above the most populated one hold 30+ regions;\na fit needs 4, so there is none",
        width=7.0,
    )
    ax = axes[0][0]
    draw_counts(ax, record, case)
    peak = max(case["bins"], key=lambda b: b[1])
    ax.annotate(
        "most populated bin:\nit and all below are left out",
        xy=peak,
        xytext=(peak[0] + 0.4, peak[1] * 1.25),
        fontsize=9,
        arrowprops={"arrowstyle": "->", "color": COUNT_COLOR},
    )
    ax.set_xticks([s for s, _ in case["bins"]])
    plain_count_ticks(ax)
    return fig


def plot_rising_test(record):
    """`test_fits_refuse_a_rising_curve`."""
    rising, against_flat, falling = record["cases"]
    handles = count_legend(reference=True)
    handles = [handles[0], handles[1], handles[2], handles[5]]
    fig, axes = figure_with_legend(
        1,
        3,
        handles,
        "Counts that rise with the score give no fit; the same counts falling do",
        width=4.4,
    )
    for ax, case in zip(axes[0], (rising, against_flat, falling)):
        draw_counts(ax, record, case)
        plain_count_ticks(ax)
        ax.set_xticks([s for s, _ in case["bins"]])
        result = fit_label(case["fit"], None) if case["fit"] else "no fit"
        ax.set_title(f"{case['label']}\n{result}", fontsize=10)
    return fig


PLOTTERS = {
    "survival_and_bin_counts": plot_counts_test,
    "fit_scores_recovers_slope_and_stops_at_homolog_excess": plot_single_curve_test,
    "fit_scores_with_reference_stops_where_the_ratio_rises": plot_reference_test,
    "fit_scores_needs_enough_bins": plot_too_few_bins_test,
    "fits_refuse_a_rising_curve": plot_rising_test,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("data_dir", help="directory the tests wrote their JSON into")
    parser.add_argument("--output", required=True, help="directory for one PNG per test")
    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    for test, plotter in PLOTTERS.items():
        with open(os.path.join(args.data_dir, f"{test}.json")) as handle:
            fig = plotter(json.load(handle))
        fig.tight_layout(rect=(0, 0, 1, fig.panels_top))
        fig.savefig(os.path.join(args.output, f"{test}.png"), dpi=130)
        plt.close(fig)


if __name__ == "__main__":
    main()
