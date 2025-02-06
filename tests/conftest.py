import os

import pytest

@pytest.fixture
def testdata():
    return os.path.join(os.path.dirname(__file__), "testdata")

@pytest.fixture
def bcl2_first25(testdata):
    return os.path.join(testdata, "bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz")


@pytest.fixture
def bcl2_all300(testdata):
    return os.path.join(testdata, "uniprotkb_BCL2_AND_model_organism_9606_2025_02_06.fasta.gz")

@pytest.fixture
def ced9(testdata):
    return os.path.join(testdata, "ced9.fasta")