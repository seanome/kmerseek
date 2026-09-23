"""The walk replay agrees with the CLI test `test_cli_search_extend_mismatch_penalty`."""

import json
import os

import pytest

from plot_xdrop_walk import plot_walk, title_line, walk_both_ways

FIXTURE = os.path.join(os.path.dirname(__file__), "testdata", "bcl2_vs_ced9.hp.k12.pair.json")


@pytest.fixture
def pair():
    with open(FIXTURE) as handle:
        return json.load(handle)


def test_bh1_seed_grows_seven_right_and_none_left(pair):
    (left_kept, left_steps), (right_kept, right_steps) = walk_both_ways(pair, pair["regions"][0])
    assert (left_kept, right_kept) == (0, 7)
    assert [s for _, _, s in left_steps] == [-2, -4, -3, -5, -7, -6, -8, -10]
    assert [s for _, _, s in right_steps] == [-2, -1, 0, -2, -1, 0, 1, -1, -3, -5, -7, -6, -8]
    # The walk starts at the seed edge and stops once it is more than 8 below its best.
    assert (left_steps[0][0], left_steps[-1][0]) == (137, 130)
    assert (right_steps[0][0], right_steps[-1][0]) == (157, 169)


def test_a_penalty_above_the_margin_stops_at_the_first_mismatch(pair):
    (left_kept, left_steps), (right_kept, right_steps) = walk_both_ways(
        pair, pair["regions"][0], penalty=9.0
    )
    assert (left_kept, right_kept) == (0, 0)
    assert ([s for _, _, s in left_steps], [s for _, _, s in right_steps]) == ([-9], [-9])


def test_title_states_both_sides(pair):
    assert title_line(pair, 0, 7, 2.0, 8.0).startswith(
        "The walk grows the seed 7 residues to the right and 0 to the left (BCL2_HUMAN vs CED9_CAEEL"
    )


def test_figure_has_the_letters_and_score_panels(pair):
    fig = plot_walk(pair, pair["regions"][0])
    assert [ax.get_ylabel() for ax in fig.axes] == ["", "running score"]
    assert fig.axes[1].get_xlabel() == "position in BCL2_HUMAN (0-based, as in the CSV)"
