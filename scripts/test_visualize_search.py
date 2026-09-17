"""Tests for visualize_search.py.

Run with: /Users/olga/anaconda3/envs/kmerseek-dev/bin/python3 -m pytest scripts/test_visualize_search.py -v

The search fixture is CED-9 (tests/testdata/fasta/ced9.fasta) searched against the
25-protein BCL-2 family FASTA at hp_lehninger2 k=12 with every filter open
(--max-query-pvalue 1 --min-region-score 0): 362 rows, 24 targets. The Pfam table is
those targets' rows from Pfam-A.regions.
"""

import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import visualize_search as vs
from visualize_hits import iter_query_rows, load_rows_for_queries, load_query_names, scan_csv

HERE = os.path.dirname(__file__)
TESTDATA = os.path.join(HERE, "testdata")
FASTA_DIR = os.path.join(HERE, "..", "tests", "testdata", "fasta")
SEARCH_CSV = os.path.join(TESTDATA, "ced9_vs_bcl2_25.hp.k12.search.csv")
CED9_FASTA = os.path.join(FASTA_DIR, "ced9.fasta")
TARGETS_FASTA = os.path.join(FASTA_DIR, "bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz")
DOMAIN_TABLES = [os.path.join(TESTDATA, "bcl2_25_pfam_domains.tsv"), os.path.join(TESTDATA, "bcl2_ced9_pfam_domains.tsv")]


@pytest.fixture
def rows():
    lazy = scan_csv(SEARCH_CSV)
    (query_name, rows), = iter_query_rows(load_rows_for_queries(lazy, load_query_names(lazy)))
    return query_name, rows


def test_read_fasta_handles_gzip_and_plain():
    ced9 = vs.read_fasta(CED9_FASTA)
    assert list(ced9) == ["sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1"]
    assert len(next(iter(ced9.values()))) == 280
    targets = vs.read_fasta(TARGETS_FASTA)
    assert len(targets) == 25
    assert len(targets["sp|P10415|BCL2_HUMAN Apoptosis regulator Bcl-2 OS=Homo sapiens OX=9606 GN=BCL2 PE=1 SV=2"]) == 239


def test_description_drops_the_token_and_the_organism():
    assert vs.description("sp|P10415|BCL2_HUMAN Apoptosis regulator Bcl-2 OS=Homo sapiens OX=9606") == "Apoptosis regulator Bcl-2"
    assert vs.description("gene1") == ""


def test_rank_targets_orders_by_q_value_then_region_score(rows):
    _, rows = rows
    ranked = vs.rank_targets(rows, max_hits=100)
    assert len(ranked) == 24
    names = [t.split("|")[2].split()[0] for t, _, _ in ranked]
    assert names[:3] == ["FBX10_HUMAN", "B2L14_HUMAN", "BAK_HUMAN"]
    q_values = [q for _, _, q in ranked]
    assert q_values == sorted(q_values)
    assert vs.rank_targets(rows, max_hits=5)[4][0] == ranked[4][0]


def test_hit_entry_reads_the_numbers_off_the_best_row(rows):
    _, rows = rows
    target, row, q = vs.rank_targets(rows, 100)[0]
    model = {"target": {"label": "FBX10_HUMAN"}, "runs": [1, 2, 3], "alphabet": "x"}
    hit = vs.hit_entry(1, target, row, q, model)
    assert hit["target"] == "FBX10_HUMAN"
    assert hit["n_runs"] == 3
    assert hit["n_shared"] == 50
    assert hit["region_score"] == pytest.approx(5.99, abs=0.01)
    assert hit["q_value"] == pytest.approx(q)


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
def test_report_end_to_end(tmp_path):
    args = vs._build_arg_parser().parse_args([
        "--csv", SEARCH_CSV, "--query-fasta", CED9_FASTA, "--target-fasta", TARGETS_FASTA,
        "--output-dir", str(tmp_path), "--domains", *DOMAIN_TABLES, "--max-hits", "12",
    ])
    lazy = scan_csv(SEARCH_CSV)
    report = vs.SearchReport(args, _kmerseek(), vs.load_domains(args.domains))
    (query_name, rows), = iter_query_rows(load_rows_for_queries(lazy, load_query_names(lazy)))
    html = report.render(query_name, rows)
    assert "<title>CED9_CAEEL kmerseek hits</title>" in html
    assert html.count('"rank":') == 12
    assert '"target": "BCL2_HUMAN"' in html
    # The query's Pfam domains label the track, and BCL-2's label its side of the pair view.
    assert '"query": {"label": "CED9_CAEEL"' in html
    assert '{"name": "BH4", "start": 76, "end": 101}' in html
    assert "Bcl-2" in html
    assert html.count("</script>") == 1
