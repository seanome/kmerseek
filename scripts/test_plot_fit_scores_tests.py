"""The figure script picks the same bins as `fit_scores_with_reference` in src/rust/evalue.rs."""

import math

import matplotlib.pyplot as plt

from plot_fit_scores_tests import PLOTTERS, plot_too_few_bins_test, ratio_baseline, usable_bins


def test_usable_bins_skip_the_peak_and_small_bins():
    bins = [[12, 100], [13, 60], [14, 40], [15, 31], [16, 29]]
    assert usable_bins(bins, [], 30) == [(13, 60, None), (14, 40, None), (15, 31, None)]


def test_usable_bins_need_the_reference_bin_too():
    bins = [[12, 100], [13, 60], [14, 40], [15, 31]]
    reference = [[12, 100], [13, 45], [14, 29]]
    assert usable_bins(bins, reference, 30) == [(13, 60, 45)]


def test_ratio_baseline_of_a_constant_ratio_is_flat():
    usable = [(s, 150 * 2**-s, 100 * 2**-s) for s in range(13, 25)]
    slope, intercept = ratio_baseline(usable, 8)
    assert abs(slope) < 1e-12
    assert abs(intercept - math.log(1.5)) < 1e-12


def test_every_test_has_a_plotter():
    assert sorted(PLOTTERS) == [
        "fit_scores_needs_enough_bins",
        "fit_scores_recovers_slope_and_stops_at_homolog_excess",
        "fit_scores_with_reference_stops_where_the_ratio_rises",
        "fits_refuse_a_rising_curve",
        "survival_and_bin_counts",
    ]


def test_too_few_bins_figure_draws_every_nonempty_bin():
    record = {
        "min_bin_count": 30,
        "cases": [
            {
                "label": "three usable bins",
                "bins": [[12, 100], [13, 60], [14, 40], [15, 31]],
                "reference_bins": [],
                "made_with_lambda": None,
                "fit": None,
            }
        ],
    }
    fig = plot_too_few_bins_test(record)
    points = [
        line.get_xydata()[0].tolist() for line in fig.axes[0].lines if line.get_marker() == "o"
    ]
    assert points == [[12, 100], [13, 60], [14, 40], [15, 31]]
    plt.close(fig)
