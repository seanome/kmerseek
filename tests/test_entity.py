import os

import polars as pl
import sourmash

from kmerseek.entity import KmerseekEntity


def test_entity(capsys, bcl2_first25, bcl2_first25_hp_k24_scaled5_sig):
    # Is there a better name for this object than "test" here?
    test = KmerseekEntity(bcl2_first25, "hp", 24, 5)
    assert test.sketch_kws == dict(moltype="hp", ksize=24, scaled=5)
    assert test.force == False
    assert test.sig == f"{bcl2_first25}.hp.k24.scaled5.sig.zip"
    assert os.path.exists(test.sig)

    test_sig = list(sourmash.load_file_as_signatures(test.sig)).pop()
    true_sig = list(
        sourmash.load_file_as_signatures(bcl2_first25_hp_k24_scaled5_sig)
    ).pop()

    assert test_sig == true_sig


def test_entity_extract_kmers(
    capsys,
    bcl2_first25,
    bcl2_first25_hp_k24_scaled5_sig,
    bcl2_first25_hp_k25_scaled5_kmers,
):
    # Is there a better name for this object than "test" here?
    test = KmerseekEntity(bcl2_first25, "hp", 24, 5, extract_kmers=True)
    assert test.sketch_kws == dict(moltype="hp", ksize=24, scaled=5)
    assert test.force == False
    assert test.sig == f"{bcl2_first25}.hp.k24.scaled5.sig.zip"
    assert os.path.exists(test.sig)

    test_sig = list(sourmash.load_file_as_signatures(test.sig)).pop()
    true_sig = list(
        sourmash.load_file_as_signatures(bcl2_first25_hp_k24_scaled5_sig)
    ).pop()
    assert test_sig == true_sig

    # Test the new kmer extraction behavior
    assert test.kmers_pq == f"{bcl2_first25}.hp.k24.scaled5.sig.zip.kmers.pq"
    assert os.path.exists(test.kmers_pq)

    test_kmers = pl.read_parquet(test.kmers_pq)
    true_kmers = pl.read_parquet(bcl2_first25_hp_k25_scaled5_kmers)

    # Sort by sequence name and kmer start position, and ignore sequence_file column since that includes
    # paths that changes
    true_kmers = true_kmers.sort(["sequence_name", "start"]).select(
        ["sequence_name", "kmer", "hashval", "encoded", "start"]
    )
    # Sort by the same columns, make sure we're comparing the same columns
    test_kmers = test_kmers.sort(["sequence_name", "start"]).select(true_kmers.columns)
    assert test_kmers.shape == (1712, 5)
    assert test_kmers.equals(true_kmers)
