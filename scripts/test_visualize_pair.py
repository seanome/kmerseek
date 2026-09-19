"""Tests for pair_model.py, visualize_pair.py and visualize_pair_html.py.

Run with: /Users/olga/anaconda3/envs/kmerseek-dev/bin/python3 -m pytest scripts/test_visualize_pair.py -v

The fixture is `kmerseek pair` output for human BCL-2 against C. elegans CED-9 at
hp_lehninger2 k=12 (tests/testdata/fasta/bcl2.fasta, ced9.fasta), plus their Pfam
domains from Pfam-A.regions. Regenerate the JSON with:

    kmerseek pair --query tests/testdata/fasta/bcl2.fasta --target tests/testdata/fasta/ced9.fasta \
        --ksize 12 --alphabet hp --output scripts/testdata/bcl2_vs_ced9.hp.k12.pair.json
"""

import copy
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import pair_model as pm
import visualize_pair as vp
from visualize_pair_html import embed_json, render_html

TESTDATA = os.path.join(os.path.dirname(__file__), "testdata")
PAIR_JSON = os.path.join(TESTDATA, "bcl2_vs_ced9.hp.k12.pair.json")
DOMAINS_TSV = os.path.join(TESTDATA, "bcl2_ced9_pfam_domains.tsv")

# BH1 in 0-based half-open coordinates: BCL-2 NWGR at 1-based 143, CED-9 SYGR at 167.
BH1 = {"query_start": 138, "query_end": 157, "target_start": 162, "target_end": 181, "length": 19}


@pytest.fixture
def pair():
    return vp.load_pair(PAIR_JSON)


@pytest.fixture
def domains():
    return pm.load_domains([DOMAINS_TSV])


@pytest.fixture
def model(pair, domains):
    return pm.build_model(pair, domains)


def test_fixture_is_bcl2_vs_ced9_at_hp_k12(pair):
    assert pair["ksize"] == 12
    assert pair["moltype"] == "hp_lehninger2"
    assert pair["query"]["name"].startswith("sp|P10415|BCL2_HUMAN")
    assert pair["target"]["name"].startswith("sp|P41958|CED9_CAEEL")
    assert len(pair["shared_kmers"]) == 27
    assert len(pair["regions"]) == 13


def test_class_residues_reads_the_lehninger_partition_off_the_sequences(pair):
    assert pm.class_residues(pair) == {"h": "AFGILMPVWY", "p": "CDEHKNQRST"}


def test_runs_drop_single_kmer_regions(pair):
    assert [r["length"] for r in pm.runs(pair)] == [19, 14, 14, 13, 13, 13]
    assert pm.runs(pair)[0] == BH1
    assert len(pm.singles(pair)) == 7
    assert (30, 94) in [(k["query_pos"], k["target_pos"]) for k in pm.singles(pair)]


def test_in_region_requires_the_same_diagonal():
    region = {"query_start": 10, "query_end": 30, "target_start": 50, "target_end": 70, "length": 20}
    assert pm.in_region({"query_pos": 10, "target_pos": 50}, region, ksize=12)
    assert pm.in_region({"query_pos": 18, "target_pos": 58}, region, ksize=12)
    assert not pm.in_region({"query_pos": 19, "target_pos": 59}, region, ksize=12)
    assert not pm.in_region({"query_pos": 12, "target_pos": 50}, region, ksize=12)


# -- domains --


def test_domain_table_reads_pfam_columns(domains):
    assert domains == [
        {"protein": "P10415", "start": 97, "end": 195, "name": "Bcl-2"},
        {"protein": "P10415", "start": 7, "end": 32, "name": "BH4"},
        {"protein": "P41958", "start": 116, "end": 221, "name": "Bcl-2"},
        {"protein": "P41958", "start": 76, "end": 101, "name": "BH4"},
    ]


def test_header_keys_cover_header_token_accession_and_entry_name():
    assert pm.header_keys("sp|P10415|BCL2_HUMAN Apoptosis regulator Bcl-2 OS=Homo sapiens") == {
        "sp|P10415|BCL2_HUMAN Apoptosis regulator Bcl-2 OS=Homo sapiens",
        "sp|P10415|BCL2_HUMAN",
        "P10415",
        "BCL2_HUMAN",
    }
    assert pm.header_keys("gene1") == {"gene1"}


def test_domains_for_sorts_by_start_and_keeps_a_named_span_over_its_accession(domains):
    rows = domains + [{"protein": "P10415", "start": 97, "end": 195, "name": "PF00452"}]
    assert pm.domains_for(rows, "sp|P10415|BCL2_HUMAN Apoptosis regulator Bcl-2") == [
        {"name": "BH4", "start": 7, "end": 32},
        {"name": "Bcl-2", "start": 97, "end": 195},
    ]
    accession_first = [{"protein": "P10415", "start": 97, "end": 195, "name": "PF00452"}] + domains
    assert [d["name"] for d in pm.domains_for(accession_first, "P10415")] == ["BH4", "Bcl-2"]


def test_region_name_takes_the_domain_with_the_most_overlap(domains):
    bcl2 = pm.domains_for(domains, "P10415")
    assert pm.region_name(bcl2, 139, 157) == "Bcl-2"
    assert pm.region_name(bcl2, 32, 45) == "BH4"  # one residue of BH4 beats no overlap
    assert pm.region_name(bcl2, 201, 213) == pm.NO_REGION
    assert pm.region_name([], 1, 10) == pm.NO_REGION


# -- model --


def test_model_labels_and_title(model):
    assert model["alphabet"] == "2-letter hydrophobic/polar alphabet (Lehninger)"
    assert [c["label"] for c in model["classes"]] == ["hydrophobic (A F G I L M P V W Y)", "polar (C D E H K N Q R S T)"]
    assert pm.title_lines(model) == [
        "BCL2_HUMAN (query) vs CED9_CAEEL (target): 27 shared 12-mers in the 2-letter hydrophobic/polar alphabet (Lehninger)",
        "a shared 12-mer is 12 consecutive residues with the same hydrophobic/polar pattern in both proteins",
    ]
    assert model["query"]["domains"] == [{"name": "BH4", "start": 7, "end": 32}, {"name": "Bcl-2", "start": 97, "end": 195}]
    assert model["target"]["length"] == 280


def test_run_blocks_are_numbered_longest_first_with_regions_and_counts(model):
    headers = [pm.run_header(b) for b in model["runs"]]
    assert headers == [
        "Run 1 · Bcl-2 × Bcl-2 · 19 aa · 5 identical · 5 of 19 polar",
        "Run 2 · BH4 × no region · 14 aa · 2 identical · 1 of 14 polar",
        "Run 3 · Bcl-2 × Bcl-2 · 14 aa · 2 identical · 6 of 14 polar",
        "Run 4 · no region × no region · 13 aa · 2 identical · 1 of 13 polar",
        "Run 5 · no region × no region · 13 aa · 1 identical · 1 of 13 polar",
        "Run 6 · no region × no region · 13 aa · 2 identical · 4 of 13 polar",
    ]
    bh1 = model["runs"][0]
    assert bh1["query_row"] == "RDGVNWGRIVAFFEFGGVM"
    assert bh1["target_row"] == "QCPMSYGRLIGLISFGGFV"
    assert bh1["middle"] == "      GR      FGG  "
    assert bh1["n_kmers"] == 8


def test_flank_adds_class_marks_to_the_middle_line(pair, domains):
    block = pm.build_model(pair, domains, flank=3)["runs"][0]
    assert block["window"] == {"query_start": 135, "query_end": 160, "target_start": 159, "target_end": 184}
    assert block["query_row"] == "ELFRDGVNWGRIVAFFEFGGVMCVE"
    # Inside the run every column is the same class, so the middle line is solid there;
    # in the flanks E/Q share a class but L/T, F/D and V/A do not.
    assert block["middle"] == ":  ::::::GR::::::FGG:: ::"
    assert pm.middle_line("AC", "AD", "hp", "hp", show_class=False) == "A "


def test_protein20_model_has_no_classes(pair):
    plain = copy.deepcopy(pair)
    plain["moltype"] = "protein20"
    for side in ("query", "target"):
        plain[side]["encoded"] = plain[side]["sequence"]
    model = pm.build_model(plain)
    assert model["classes"] == []
    assert model["alphabet"] == "protein20 (no reduction)"
    assert model["runs"][0]["polar"] is None
    assert pm.run_header(model["runs"][0]) == "Run 1 · no region × no region · 19 aa · 5 identical"
    assert pm.title_lines(model)[1] == "a shared 12-mer is 12 consecutive identical residues in both proteins"


# -- renderers --


def test_output_basename(pair):
    assert vp.output_basename(pair) == "BCL2_HUMAN_vs_CED9_CAEEL.hp_lehninger2.k12"


def test_figure_writes_png_and_svg_with_every_run_and_legend_entry(model, tmp_path):
    png, svg = tmp_path / "pair.png", tmp_path / "pair.svg"
    vp.plot_pair(model, [str(png), str(svg)])
    assert png.stat().st_size > 0
    text = svg.read_text()
    for expected in (
        "hydrophobic (A F G I L M P V W Y)",
        "run of 2 or more consecutive shared 12-mers (6), numbered",
        "single shared 12-mer (7)",
        "protein, with its domains as boxes; each domain's span shaded across the plot",
        "Run 1 · Bcl-2 × Bcl-2 · 19 aa · 5 identical · 5 of 19 polar",
        "Run 6 · no region × no region · 13 aa · 2 identical · 4 of 13 polar",
        "BCL2_HUMAN position (aa)",
        "CED9_CAEEL position (aa)",
        "BH4",
    ):
        assert expected in text, expected


def test_figure_without_runs_or_domains(pair, tmp_path):
    singles_only = copy.deepcopy(pair)
    singles_only["regions"] = [r for r in pair["regions"] if r["length"] == 12]
    model = pm.build_model(singles_only)
    assert model["runs"] == [] and len(model["singles"]) == 27
    svg = tmp_path / "pair.svg"
    vp.plot_pair(model, [str(svg)])
    text = svg.read_text()
    assert "single shared 12-mer (27)" in text
    assert "protein, with its domains" not in text


def test_html_embeds_the_model_once(model, tmp_path):
    path = tmp_path / "pair.html"
    vp.write_html(model, str(path))
    text = path.read_text()
    assert "<title>BCL2_HUMAN vs CED9_CAEEL shared k-mers</title>" in text
    assert '"query_row": "RDGVNWGRIVAFFEFGGVM"' in text
    assert text.count("<script>") == 1 and text.count("</script>") == 1


def test_embed_json_never_lets_a_name_close_the_script():
    out = embed_json({"name": "x</script><b>"})
    assert "</script>" not in out
    assert "x<\\/script>" in out
    page = render_html({"query": {"label": "x</script>"}, "target": {"label": "y"}})
    assert page.count("</script>") == 1


# -- gapped identities and the structural path --

USALIGN_REPORT = os.path.join(TESTDATA, "bcl2_vs_ced9.usalign.txt")


@pytest.fixture
def structure():
    import structure_alignment as sa

    with open(USALIGN_REPORT) as fh:
        report = sa.parse_report(fh.read())
    return report | {"aligner": "USalign", "query_file": "AF-P10415-F1-model_v6.cif", "target_file": "AF-P41958-F1-model_v6.cif"}


def test_gapped_block_counts_identities_over_the_run_columns(pair, domains):
    model = pm.build_model(pair, domains, gap_flank=10)
    bh1 = model["runs"][0]["gapped"]
    # Run 1 with ten residues either side aligns without a gap; the run's 19 columns hold
    # the same 5 identities as its diagonal.
    assert (bh1["query_start"], bh1["query_end"], bh1["target_start"], bh1["target_end"]) == (128, 167, 152, 191)
    assert bh1["run_columns"] == [10, 29]
    assert bh1["query_row"][10:29] == "RDGVNWGRIVAFFEFGGVM"
    assert (bh1["identical"], bh1["aligned"]) == (5, 19)
    assert bh1["middle"][10:29] == "::::::GR::::::FGG::"
    assert pm.run_header(model["runs"][0]).startswith(
        "Run 1 · Bcl-2 × Bcl-2 · 19 aa · 5 identical on the run's diagonal, 5 of its 19 aligned columns after gapped alignment"
    )
    # Run 6 sits at CED-9's C terminus, so the window is short there and the alignment gaps.
    run6 = model["runs"][5]["gapped"]
    assert "-" in run6["query_row"] + run6["target_row"]
    assert run6["target_enc"].count("-") == run6["target_row"].count("-")


def test_structure_offset_puts_run_1_on_the_path_and_run_3_off_it(pair, domains, structure):
    model = pm.build_model(pair, domains, structure=structure)
    offsets = [b["structure_offset"] for b in model["runs"]]
    assert offsets == [0, -157, 4, -166, -154, -40]
    assert [pm.structure_phrase(o) for o in offsets[:3]] == [
        "on the structural path",
        "157 residues off the structural path",
        "4 residues off the structural path",
    ]
    assert pm.structure_phrase(None) == "not structurally aligned"
    assert pm.title_lines(model)[2] == (
        "USalign of AF-P10415-F1-model_v6.cif against AF-P41958-F1-model_v6.cif: TM-score 0.55 (by BCL2_HUMAN length), "
        "RMSD 2.8 Å over 154 aligned residues"
    )
    assert "structure_offset" not in pm.build_model(pair, domains)["runs"][0]


def test_figure_and_html_draw_the_structural_path(pair, domains, structure, tmp_path):
    model = pm.build_model(pair, domains, gap_flank=10, structure=structure)
    svg = tmp_path / "pair.svg"
    vp.plot_pair(model, [str(svg)])
    text = svg.read_text()
    assert "structural alignment (USalign, TM-score 0.55): every aligned residue pair" in text
    assert "on the structural path" in text and "4 residues off the structural path" in text
    assert vp.path_segments([(0, 0, True), (1, 1, True), (5, 9, False)]) == [([0, 1], [0, 1]), ([5], [9])]
    page = render_html(model)
    assert '"tm_score_query": 0.55352' in page and '"structure_offset": 4' in page
