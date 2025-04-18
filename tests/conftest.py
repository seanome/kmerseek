import os

import pytest


@pytest.fixture
def testdata():
    return os.path.join(os.path.dirname(__file__), "testdata")


# --- Fasta files --- #
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
def bcl2_first25_hp_k24_scaled5_sig(bcl2_first25):
    return f"{bcl2_first25}.hp.k24.scaled5.sig.zip"


@pytest.fixture
def bcl2_first25_hp_k25_scaled5_kmers(bcl2_first25):
    return f"{bcl2_first25}.hp.k24.scaled5.sig.TRUE.zip.kmers.pq"


@pytest.fixture
def bcl2_all300(fasta_folder):
    return os.path.join(
        fasta_folder, "uniprotkb_BCL2_AND_model_organism_9606_2025_02_06.fasta.gz"
    )


@pytest.fixture
def ced9(fasta_folder):
    return os.path.join(fasta_folder, "ced9.fasta")


# --- Indexed Rocksdb files --- #
@pytest.fixture
def index_folder(testdata):
    return os.path.join(testdata, "index")


@pytest.fixture
def bcl2_rocksdb(index_folder):
    return os.path.join(
        index_folder,
        "bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz.hp.k15.scaled5.sig.zip.rocksdb",
    )


@pytest.fixture
def bcl2_hp_k16_sig_zip(index_folder):
    return os.path.join(
        index_folder,
        "bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz.hp.k16.scaled5.sig.zip",
    )
