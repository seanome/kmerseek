"""Tests for gapped_alignment.py.

Run with: /Users/olga/anaconda3/envs/kmerseek-dev/bin/python3 -m pytest scripts/test_gapped_alignment.py -v
"""

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from gapped_alignment import BLOSUM62, GAP_EXTEND, GAP_OPEN, global_align, local_align


def test_blosum62_is_symmetric_with_the_known_diagonal():
    assert BLOSUM62.shape == (24, 24)
    assert (BLOSUM62 == BLOSUM62.T).all()
    # W/W 11, C/C 9, A/A 4, W/C -2, the values everyone checks first.
    assert BLOSUM62[17, 17] == 11 and BLOSUM62[4, 4] == 9 and BLOSUM62[0, 0] == 4 and BLOSUM62[17, 4] == -2


def test_identical_sequences_align_without_gaps():
    r = global_align("NWGRIVAFFEFGG", "NWGRIVAFFEFGG")
    assert r["a_row"] == r["b_row"] == "NWGRIVAFFEFGG"
    assert (r["identical"], r["aligned"], r["columns"]) == (13, 13, 13)
    assert r["score"] == sum(int(BLOSUM62[i, i]) for i in map("ARNDCQEGHILKMFPSTWYV".index, "NWGRIVAFFEFGG"))


def test_global_alignment_charges_end_gaps():
    r = global_align("AAAA", "WW")
    assert r["a_row"] == "AAAA" and r["b_row"] == "--WW"
    assert r["score"] == -(GAP_OPEN + GAP_EXTEND) + 2 * int(BLOSUM62[0, 17])
    assert (r["identical"], r["aligned"], r["columns"]) == (0, 2, 4)


def test_bcl2_against_mcl1_bh1_needs_one_gap_to_line_up_nwgr():
    # BCL-2 129-153 (1-based) and MCL-1 245-273: the exact hp run sits 3 residues off the
    # true alignment, which one gap restores.
    bcl2 = "ELFRDGVNWGRIVAFFEFGGVMCVE"
    mcl1 = "VMIHVFSDGVTNWGRIVTLISFGAFVAKH"
    r = global_align(bcl2, mcl1)
    assert r["a_row"] == "---ELFRDGV-NWGRIVAFFEFGGVMCVE"
    assert r["b_row"] == "VMIHVFSDGVTNWGRIVTLISFGAFVAKH"
    assert r["identical"] == 12
    assert r["aligned"] == 25


def test_local_alignment_trims_to_the_scoring_core():
    r = local_align("HEAGAWGHEE", "PAWHEAE")
    assert (r["a_row"], r["b_row"], r["score"]) == ("HEA", "HEA", 17)
    assert (r["a_start"], r["a_end"], r["b_start"], r["b_end"]) == (0, 3, 3, 6)
    assert local_align("AAAA", "WWWW") is None
    assert local_align("", "WWWW") is None
