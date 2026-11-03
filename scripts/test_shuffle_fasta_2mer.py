"""The 2-mer shuffle keeps every dipeptide count, and the Rust decoy fixture regenerates from its seed."""

import gzip
import os
import random
import subprocess
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from shuffle_fasta_2mer import parse_fasta, shuffle_2mer_python  # noqa: E402

FASTA = os.path.join(ROOT, "tests/testdata/fasta")
QUERIES = os.path.join(FASTA, "bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz")
DECOYS = os.path.join(FASTA, "bcl2_25_shuffled_2mer_20x_seed1.fasta.gz")
# Human BCL2 residues 1-30 (UniProt P10415), the BCL2_FRAGMENT_1_30 of the Rust tests.
BCL2_1_30 = "MAHAGRTGYDNREIVMKYIHYKLSQRGYEW"


def dipeptides(seq):
    return Counter(zip(seq, seq[1:]))


def test_every_dipeptide_count_and_both_ends_are_kept():
    rng = random.Random(1)
    for _ in range(200):
        shuffled = shuffle_2mer_python(BCL2_1_30, rng)
        assert dipeptides(shuffled) == dipeptides(BCL2_1_30)
        assert (shuffled[0], shuffled[-1]) == (BCL2_1_30[0], BCL2_1_30[-1])


def test_shuffles_differ_from_the_original():
    rng = random.Random(1)
    shuffles = {shuffle_2mer_python(BCL2_1_30, rng) for _ in range(200)}
    assert len(shuffles) > 100


def test_decoy_fixture_regenerates_from_seed_1():
    out = subprocess.run(
        [sys.executable, os.path.join(ROOT, "shuffle_fasta_2mer.py"), QUERIES, "-n", "20", "--seed", "1"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    with gzip.open(DECOYS, "rt") as handle:
        assert out == handle.read()


def test_every_decoy_keeps_its_protein_dipeptides():
    proteins = dict(parse_fasta(QUERIES))
    decoys = list(parse_fasta(DECOYS))
    assert len(proteins) == 25 and len(decoys) == 500
    for header, seq in decoys:
        original = proteins[header.rsplit("_shuffle", 1)[0]]
        assert dipeptides(seq) == dipeptides(original)
