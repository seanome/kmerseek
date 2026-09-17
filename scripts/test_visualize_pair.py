"""Tests for visualize_pair.py.

Run with: /Users/olga/anaconda3/envs/kmerseek-dev/bin/python3 -m pytest scripts/test_visualize_pair.py -v

The fixture is `kmerseek pair` output for human BCL-2 against C. elegans CED-9 at
hp_lehninger2 k=12 (tests/testdata/fasta/bcl2.fasta, ced9.fasta). Regenerate it with:

    kmerseek pair --query tests/testdata/fasta/bcl2.fasta --target tests/testdata/fasta/ced9.fasta \
        --ksize 12 --alphabet hp --output scripts/testdata/bcl2_vs_ced9.hp.k12.pair.json
"""

import copy
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import visualize_pair as vp

PAIR_JSON = os.path.join(os.path.dirname(__file__), "testdata", "bcl2_vs_ced9.hp.k12.pair.json")

# BH1 in 0-based half-open coordinates: BCL-2 NWGR at 1-based 143, CED-9 SYGR at 167.
BH1 = {"query_start": 138, "query_end": 157, "target_start": 162, "target_end": 181, "length": 19}


@pytest.fixture
def pair():
    return vp.load_pair(PAIR_JSON)


def test_fixture_is_bcl2_vs_ced9_at_hp_k12(pair):
    assert pair["ksize"] == 12
    assert pair["moltype"] == "hp_lehninger2"
    assert pair["query"]["name"].startswith("sp|P10415|BCL2_HUMAN")
    assert pair["target"]["name"].startswith("sp|P41958|CED9_CAEEL")
    assert len(pair["shared_kmers"]) == 27
    assert len(pair["regions"]) == 14


def test_class_residues_reads_the_lehninger_partition_off_the_sequences(pair):
    assert vp.class_residues(pair) == {"h": "AFGILMPVWY", "p": "CDEHKNQRST"}


def test_runs_drop_single_kmer_regions(pair):
    run_lengths = [r["length"] for r in vp.runs(pair)]
    assert run_lengths == [19, 14, 14, 13, 13]
    assert vp.longest_run(pair) == BH1


def test_ribbon_window_adds_flank_on_both_sides(pair):
    window = vp.ribbon_window(BH1, 239, 280, flank=10)
    assert window == {"query_start": 128, "query_end": 167, "target_start": 152, "target_end": 191}
    slices = vp.window_slices(pair, window)
    assert slices["query_seq"] == "RFATVVEELFRDGVNWGRIVAFFEFGGVMCVESVNREMS"
    assert slices["target_seq"] == "VRTVGNAQTDQCPMSYGRLIGLISFGGFVAAKMMESVEL"
    assert slices["query_enc"] == "phhphhpphhpphhphhphhhhhphhhhhphpphppphp"


def test_ribbon_window_clips_at_the_sequence_ends():
    region = {"query_start": 2, "query_end": 20, "target_start": 5, "target_end": 23, "length": 18}
    window = vp.ribbon_window(region, query_len=25, target_len=100, flank=10)
    assert window == {"query_start": 0, "query_end": 25, "target_start": 3, "target_end": 28}


def test_region_counts_match_the_figure(pair):
    assert vp.region_counts(pair, BH1) == {"identical": 5, "same_class": 19, "length": 19}


def test_split_kmers_by_run(pair):
    on_run, single = vp.split_kmers_by_run(pair)
    assert len(on_run) == 18
    assert len(single) == 9
    bh1_kmers = [(k["query_pos"], k["target_pos"]) for k in on_run if 138 <= k["query_pos"] < 157]
    assert bh1_kmers == [(138 + i, 162 + i) for i in range(8)]
    assert (30, 94) in [(k["query_pos"], k["target_pos"]) for k in single]


def test_in_region_requires_the_same_diagonal():
    region = {"query_start": 10, "query_end": 30, "target_start": 50, "target_end": 70, "length": 20}
    assert vp.in_region({"query_pos": 10, "target_pos": 50}, region, ksize=12)
    assert vp.in_region({"query_pos": 18, "target_pos": 58}, region, ksize=12)
    assert not vp.in_region({"query_pos": 19, "target_pos": 59}, region, ksize=12)
    assert not vp.in_region({"query_pos": 12, "target_pos": 50}, region, ksize=12)


def test_row_labels_are_one_based_inclusive(pair):
    plot = vp.PairPlot(pair)
    assert plot.row_labels(plot.window()) == [
        "CED9_CAEEL 153–191",
        "hp_lehninger2",
        "hp_lehninger2",
        "BCL2_HUMAN 129–167",
    ]


def test_title_lines(pair):
    assert vp.PairPlot(pair).title_lines() == [
        "BCL2_HUMAN (query) vs CED9_CAEEL (target), hp_lehninger2 12-mers",
        "27 shared: 18 in runs of consecutive k-mers, 9 singles; longest run 19 residues, "
        "5/19 identical, 19/19 same class",
    ]


def test_output_basename(pair):
    assert vp.output_basename(pair) == "BCL2_HUMAN_vs_CED9_CAEEL.hp_lehninger2.k12"


def test_plot_writes_png_and_svg_with_legend_text(pair, tmp_path):
    png, svg = tmp_path / "pair.png", tmp_path / "pair.svg"
    vp.plot_pair(pair, [str(png), str(svg)])
    assert png.stat().st_size > 0
    text = svg.read_text()
    for expected in (
        "class h: A F G I L M P V W Y",
        "same class in both",
        "8 consecutive shared 12-mers: BCL2_HUMAN 139–157, CED9_CAEEL 163–181",
        "single shared 12-mer (9)",
        "shared 12-mer in a run (18)",
        "run of consecutive shared k-mers (5)",
        "BCL2_HUMAN position (query, 239 aa)",
        "CED9_CAEEL position (target, 280 aa)",
    ):
        assert expected in text, expected


def test_plot_without_a_run_still_draws_the_dot_plot(pair, tmp_path):
    singles_only = copy.deepcopy(pair)
    singles_only["regions"] = [r for r in pair["regions"] if r["length"] == 12]
    plot = vp.PairPlot(singles_only)
    assert plot.region is None
    assert plot.title_lines()[1] == "27 shared, none consecutive in both sequences"
    svg = tmp_path / "pair.svg"
    vp.plot_pair(singles_only, [str(svg)])
    text = svg.read_text()
    assert "No run to show: no two shared 12-mers are consecutive in both sequences." in text
    assert "single shared 12-mer (27)" in text


def test_protein20_pair_has_no_encoded_rows(pair):
    plain = copy.deepcopy(pair)
    plain["moltype"] = "protein20"
    for side in ("query", "target"):
        plain[side]["encoded"] = plain[side]["sequence"]
    plot = vp.PairPlot(plain)
    assert not plot.reduced
    assert plot.classes == {}
    assert plot.n_rows() == 2
    assert plot.row_labels(plot.window()) == ["CED9_CAEEL 153–191", "BCL2_HUMAN 129–167"]


def test_html_embeds_the_pair_and_escapes_script_closers(pair, tmp_path):
    path = tmp_path / "pair.html"
    vp.write_html(pair, str(path), flank=7)
    text = path.read_text()
    assert "<title>BCL2_HUMAN vs CED9_CAEEL shared k-mers</title>" in text
    assert 'value="7"' in text
    assert "const state = { run: 0, flank: 7 };" in text
    assert '"query_pos": 138, "target_pos": 162, "kmer": "pphhphhphhhh"' in text
    assert text.count("<script>") == 1 and text.count("</script>") == 1


def test_html_never_lets_a_sequence_name_close_the_script():
    from visualize_pair_html import render_html

    hostile = {
        "ksize": 3,
        "moltype": "protein20",
        "query": {"name": "x</script><b>", "sequence": "APG", "encoded": "APG"},
        "target": {"name": "y", "sequence": "APG", "encoded": "APG"},
        "shared_kmers": [],
        "regions": [],
    }
    text = render_html(hostile, "t")
    assert "x</script>" not in text
    assert "x<\\/script>" in text
    assert text.count("</script>") == 1
