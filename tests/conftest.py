import os

import pytest


@pytest.fixture
def testdata():
    return os.path.join(os.path.dirname(__file__), "testdata")


@pytest.fixture
def fasta_folder(testdata):
    return os.path.join(testdata, "fasta")


@pytest.fixture
def bcl2_first25(fasta_folder):
    return os.path.join(
        fasta_folder,
        "bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz",
    )


@pytest.fixture
def bcl2_all300(fasta_folder):
    return os.path.join(
        fasta_folder, "uniprotkb_BCL2_AND_model_organism_9606_2025_02_06.fasta.gz"
    )


@pytest.fixture
def ced9(fasta_folder):
    return os.path.join(fasta_folder, "ced9.fasta")
