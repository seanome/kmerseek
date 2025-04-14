import os
import polars as pl
from click.testing import CliRunner

from kmerseek.main import cli


def test_index(bcl2_first25):
    runner = CliRunner()
    result = runner.invoke(cli, ["index", "--force", bcl2_first25])
    assert result.exit_code == 0

    # Make sure all the files got created
    assert os.path.exists(f"{bcl2_first25}.manysketch.csv")
    with open(f"{bcl2_first25}.manysketch.csv") as f:
        assert f.readlines() == [
            "name,genome_filename,protein_filename\n",
            f"{os.path.basename(bcl2_first25)},,{bcl2_first25}\n",
        ]

    # TODO: Actually test for the signature in this file by reading it in with
    # `sourmash.load_file_as_signatures`
    # https://github.com/seanome/kmerseek/issues/3
    assert os.path.exists(f"{bcl2_first25}.hp.k24.scaled5.sig.zip")

    assert os.path.exists(f"{bcl2_first25}.hp.k24.scaled5.sig.zip.siglist")
    with open(f"{bcl2_first25}.hp.k24.scaled5.sig.zip.siglist") as f:
        assert f.readlines() == [f"{bcl2_first25}.hp.k24.scaled5.sig.zip"]

    # TODO: Actually test this, e.g. with what Sourmash Branchwater does for testing indices
    # https://github.com/sourmash-bio/sourmash_plugin_branchwater/blob/776727091d78af7de4990578ec6d07b12882c1b8/src/python/tests/test_index.py#L595
    # https://github.com/seanome/kmerseek/issues/3
    assert os.path.exists(f"{bcl2_first25}.hp.k24.scaled5.sig.zip.rocksdb")


def test_index_extract_kmers(bcl2_first25):
    runner = CliRunner()
    result = runner.invoke(
        cli, ["index", "--extract-kmers", "--force", bcl2_first25]
    )
    assert result.exit_code == 0

    # Make sure all the files got created
    assert os.path.exists(f"{bcl2_first25}.manysketch.csv")
    with open(f"{bcl2_first25}.manysketch.csv") as f:
        assert f.readlines() == [
            "name,genome_filename,protein_filename\n",
            f"{os.path.basename(bcl2_first25)},,{bcl2_first25}\n",
        ]

    # TODO: Actually test for the signature in this file by reading it in with
    # `sourmash.load_file_as_signatures`
    # https://github.com/seanome/kmerseek/issues/3
    assert os.path.exists(f"{bcl2_first25}.hp.k24.scaled5.sig.zip")

    # TODO: Actually test this by reading it in with polars, then:
    # - assert the `.head()`
    # - assert the `.tail()`
    # - assert the `.describe()`
    # May need to `.sort_values()` to make sure `head()` and `tail()` work
    # https://github.com/seanome/kmerseek/issues/3
    assert os.path.exists(f"{bcl2_first25}.hp.k24.scaled5.sig.zip.kmers.pq")
    kmers_true = pl.read_parquet(f"{bcl2_first25}.hp.k24.scaled5.sig.TRUE.zip.kmers.pq")
    # Sort by sequence name and kmer start position
    kmers_true = kmers_true.sort(["sequence_name", "start"]).select(
        ["sequence_name", "kmer", "hashval", "encoded", "start"]
    )

    kmers_test = pl.read_parquet(f"{bcl2_first25}.hp.k24.scaled5.sig.zip.kmers.pq")
    # Sort by sequence name and kmer start position, then make sure columns are in the same order
    kmers_test = kmers_test.sort(["sequence_name", "start"]).select(kmers_true.columns)
    assert kmers_test.shape == (1712, 5)
    assert kmers_test.equals(kmers_true)

    assert os.path.exists(f"{bcl2_first25}.hp.k24.scaled5.sig.zip.siglist")
    with open(f"{bcl2_first25}.hp.k24.scaled5.sig.zip.siglist") as f:
        assert f.readlines() == [f"{bcl2_first25}.hp.k24.scaled5.sig.zip"]

    # TODO: Actually test this, e.g. with what Sourmash Branchwater does for testing indices
    # https://github.com/sourmash-bio/sourmash_plugin_branchwater/blob/776727091d78af7de4990578ec6d07b12882c1b8/src/python/tests/test_index.py#L595
    # https://github.com/seanome/kmerseek/issues/3
    assert os.path.exists(f"{bcl2_first25}.hp.k24.scaled5.sig.zip.rocksdb")
