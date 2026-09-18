"""Tests for visualize_search.py.

Run with: /Users/olga/anaconda3/envs/kmerseek-dev/bin/python3 -m pytest scripts/test_visualize_search.py -v

The search fixture is CED-9 (tests/testdata/fasta/ced9.fasta) searched against the
25-protein BCL-2 family FASTA at hp_lehninger2 k=12 with every filter open
(--max-query-pvalue 1 --min-region-score 0): 362 rows, 24 targets. The Pfam table is
those targets' rows from Pfam-A.regions.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import visualize_search as vs
from visualize_hits import iter_query_rows, load_query_names, load_rows_for_queries, scan_csv

HERE = os.path.dirname(__file__)
TESTDATA = os.path.join(HERE, "testdata")
FASTA_DIR = os.path.join(HERE, "..", "tests", "testdata", "fasta")
SEARCH_CSV = os.path.join(TESTDATA, "ced9_vs_bcl2_25.hp.k12.search.csv")
CED9_FASTA = os.path.join(FASTA_DIR, "ced9.fasta")
TARGETS_FASTA = os.path.join(FASTA_DIR, "bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz")
DOMAIN_TABLES = [os.path.join(TESTDATA, "bcl2_25_pfam_domains.tsv"), os.path.join(TESTDATA, "bcl2_ced9_pfam_domains.tsv")]

BCL2 = "sp|P10415|BCL2_HUMAN Apoptosis regulator Bcl-2 OS=Homo sapiens OX=9606 GN=BCL2 PE=1 SV=2"
BCL2_TREMBL = "tr|A9QXG9|A9QXG9_HUMAN Apoptosis regulator Bcl-2 OS=Homo sapiens OX=9606 GN=bcl-2 PE=2 SV=1"


@pytest.fixture
def rows():
    lazy = scan_csv(SEARCH_CSV)
    ((query_name, rows),) = iter_query_rows(load_rows_for_queries(lazy, load_query_names(lazy)))
    return query_name, rows


def test_read_fasta_handles_gzip_and_plain():
    ced9 = vs.read_fasta(CED9_FASTA)
    assert list(ced9) == ["sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1"]
    assert len(next(iter(ced9.values()))) == 280
    targets = vs.read_fasta(TARGETS_FASTA)
    assert len(targets) == 25
    assert len(targets[BCL2]) == 239


def test_header_fields_and_description():
    assert vs.header_field(BCL2, "GN") == "BCL2"
    assert vs.header_field(BCL2, "OS") == "Homo sapiens"
    assert vs.header_field("gene1", "GN") is None
    assert vs.description(BCL2) == "Apoptosis regulator Bcl-2"
    assert vs.description("gene1") == ""


def test_protein_key_folds_trembl_copies_of_one_gene():
    assert vs.protein_key(BCL2) == ("bcl2", "Homo sapiens")
    assert vs.protein_key(BCL2_TREMBL) == vs.protein_key(BCL2)
    assert vs.protein_key("gene1 some description") == ("gene1 some description", None)


def test_ranking_falls_back_to_q_value_without_region_evalue(rows):
    _, rows = rows
    stat, is_evalue = vs.ranking(rows)
    assert not is_evalue
    assert len(stat) == 24
    fbx10 = next(t for t in stat if "FBX10_HUMAN" in t)
    assert stat[fbx10] == pytest.approx(1.36e-5, rel=0.05)


def test_rank_proteins_orders_by_statistic_and_caps(rows):
    _, rows = rows
    ranked = vs.rank_proteins(rows, max_rows=100)
    assert len(ranked) == 24  # every fixture target is a different gene
    names = [t.split("|")[2].split()[0] for t, _, _ in ranked]
    assert names[:3] == ["FBX10_HUMAN", "B2L14_HUMAN", "BAK_HUMAN"]
    assert [others for _, others, _ in ranked] == [[]] * 24
    stats = [s for _, _, s in ranked]
    assert stats == sorted(stats)
    assert vs.rank_proteins(rows, max_rows=5)[4][0] == ranked[4][0]


def test_coverage_counts_entries_per_query_residue(rows):
    _, rows = rows
    cov = vs.coverage(rows, query_length=280, ksize=12, solid_identical=5)
    assert len(cov["any"]) == 280 and len(cov["solid"]) == 280
    assert max(cov["any"]) == 12
    assert cov["any"].index(max(cov["any"])) == 173  # 0-based: residue 174, inside BH1
    assert max(cov["solid"]) == 3
    assert all(s <= a for s, a in zip(cov["solid"], cov["any"]))
    assert cov["any"][:3] == [2, 3, 5]


def test_has_pair_rejects_a_binary_without_the_subcommand(tmp_path):
    fake = tmp_path / "kmerseek"
    fake.write_text("#!/bin/sh\nexit 1\n")
    fake.chmod(0o755)
    assert not vs.has_pair(str(fake))
    assert not vs.has_pair(str(tmp_path / "missing"))


def _kmerseek():
    try:
        return vs.find_kmerseek(None)
    except SystemExit:
        return None


@pytest.mark.skipif(_kmerseek() is None, reason="needs a kmerseek build with `pair`")
def test_report_end_to_end(rows, tmp_path):
    query_name, rows = rows
    args = vs._build_arg_parser().parse_args([
        "--csv", SEARCH_CSV, "--query-fasta", CED9_FASTA, "--target-fasta", TARGETS_FASTA,
        "--output-dir", str(tmp_path), "--domains", *DOMAIN_TABLES, "--max-rows", "12",
    ])
    report = vs.SearchReport(args, _kmerseek(), vs.load_domains(args.domains)).report(query_name, rows)
    assert report["query"]["label"] == "CED9_CAEEL"
    assert report["query"]["domains"] == [{"name": "BH4", "start": 76, "end": 101}, {"name": "Bcl-2", "start": 116, "end": 221}]
    assert report["stat_name"] == "q-value"
    assert (report["n_entries"], report["n_proteins"], len(report["rows"])) == (24, 24, 12)
    bcl2 = next(r for r in report["rows"] if r["label"] == "BCL2_HUMAN")
    assert (bcl2["length"], bcl2["n_runs"], bcl2["best_length"], bcl2["best_identical"], bcl2["n_shared"]) == (239, 6, 19, 5, 24)
    assert bcl2["runs"][0] == {"number": 1, "query_start": 162, "query_end": 181, "target_start": 138, "target_end": 157, "length": 19, "identical": 5, "polar": 5}
    html = vs.render_report(report)
    assert "<title>CED9_CAEEL kmerseek hits</title>" in html
    assert html.count("</script>") == 1
