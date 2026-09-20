"""Tests for gapped_alignment.py.

Run with: /Users/olga/anaconda3/envs/kmerseek-dev/bin/python3 -m pytest scripts/test_gapped_alignment.py -v
"""

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from gapped_alignment import global_align


def test_identical_sequences_align_without_gaps():
    r = global_align("NWGRIVAFFEFGG", "NWGRIVAFFEFGG")
    assert r["a_row"] == r["b_row"] == "NWGRIVAFFEFGG"
    assert (r["identical"], r["aligned"], r["columns"], r["score"]) == (13, 13, 13, 75)
    assert (r["a_start"], r["a_end"], r["b_start"], r["b_end"]) == (0, 13, 0, 13)


def test_global_alignment_charges_end_gaps():
    # One gap of two (open 11 + extend 1) and two A/W columns at -3 each.
    r = global_align("AAAA", "WW")
    assert r["a_row"] == "AAAA" and r["b_row"] == "--WW"
    assert r["score"] == -12 - 6
    assert (r["identical"], r["aligned"], r["columns"]) == (0, 2, 4)


def test_bcl2_against_mcl1_bh1_needs_one_gap_to_line_up_nwgr():
    # BCL-2 129-153 (1-based) and MCL-1 245-273: the exact hp run sits 3 residues off the
    # true alignment, which one gap restores.
    bcl2 = "ELFRDGVNWGRIVAFFEFGGVMCVE"
    mcl1 = "VMIHVFSDGVTNWGRIVTLISFGAFVAKH"
    r = global_align(bcl2, mcl1)
    assert r["a_row"] == "---ELFRDGV-NWGRIVAFFEFGGVMCVE"
    assert r["b_row"] == "VMIHVFSDGVTNWGRIVTLISFGAFVAKH"
    assert (r["identical"], r["aligned"], r["score"]) == (12, 25, 44)


def test_unknown_residues_and_lower_case_align():
    # Scored as X against X (-1) but shown as written.
    r = global_align("nwgru", "NWGRO")
    assert r["a_row"] == "nwgru" and r["b_row"] == "NWGRO"
    assert (r["identical"], r["score"]) == (0, 6 + 11 + 6 + 5 - 1)
    assert global_align("", "WWWW") is None
