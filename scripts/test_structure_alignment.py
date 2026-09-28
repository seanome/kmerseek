"""Tests for structure_alignment.py.

Run with: /Users/olga/anaconda3/envs/kmerseek-dev/bin/python3 -m pytest scripts/test_structure_alignment.py -v

The fixture is USalign's report for the AlphaFold models of BCL-2 (P10415) and CED-9
(P41958), `USalign AF-P10415-F1-model_v6.cif AF-P41958-F1-model_v6.cif -outfmt 0`.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import structure_alignment as sa

REPORT = os.path.join(os.path.dirname(__file__), "testdata", "bcl2_vs_ced9.usalign.txt")
BCL2 = "sp|P10415|BCL2_HUMAN Apoptosis regulator Bcl-2 OS=Homo sapiens"


def test_parse_report_reads_scores_and_pairs():
    with open(REPORT) as fh:
        r = sa.parse_report(fh.read())
    assert (r["tm_score_query"], r["tm_score_target"], r["rmsd"], r["aligned"], r["seq_id"]) == (0.55352, 0.47977, 2.8, 154, 0.188)
    assert len(r["pairs"]) == 154
    partner = {q: (t, close) for q, t, close in r["pairs"]}
    # BH1: BCL-2 NWGR (0-based 142..145) pairs with CED-9 SYGR (166..169), within 5 A.
    assert partner[142] == (166, True) and partner[145] == (169, True)
    # BCL-2's first aligned residue is its third, paired far from CED-9's 68th.
    assert min(partner) == 2 and partner[2] == (67, False)
    # The loop (BCL-2 residues 40..85, 0-based) has no structural partner.
    assert all(q not in partner for q in range(40, 85))


def test_residue_pairs_walk_gaps_in_both_rows():
    # Column 0 gaps b, column 2 gaps a; marks are ':' close, '.' far, ' ' unaligned.
    pairs = sa.residue_pairs("MA-KL", "   :.", "-AWKL")
    assert pairs == [(1, 0, False), (2, 2, True), (3, 3, False)]


def test_find_structure_by_accession_prefers_the_newest_model(tmp_path):
    for name in ("AF-P10415-F1-model_v4.pdb", "AF-P10415-F1-model_v6.cif", "AF-Q07817-F1-model_v6.cif"):
        (tmp_path / name).write_text("")
    assert sa.find_structure(str(tmp_path), BCL2).endswith("AF-P10415-F1-model_v6.cif")
    assert sa.find_structure(str(tmp_path), "sp|Q07817|B2CL1_HUMAN").endswith("AF-Q07817-F1-model_v6.cif")
    assert sa.find_structure(str(tmp_path), "sp|P41958|CED9_CAEEL") is None
    (tmp_path / "gene1.pdb").write_text("")
    assert sa.find_structure(str(tmp_path), "gene1 some description").endswith("gene1.pdb")


def test_align_pair_returns_none_without_a_structure(tmp_path):
    assert sa.align_pair("/usr/bin/true", str(tmp_path), BCL2, "sp|P41958|CED9_CAEEL") is None
    assert sa.align_pair(None, str(tmp_path), BCL2, "sp|P41958|CED9_CAEEL") is None
