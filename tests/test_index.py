import os

from click.testing import CliRunner
import polars as pl

from kmerseek.main import cli


def test_index(bcl2_first25):
    runner = CliRunner()
    result = runner.invoke(cli, ["index", bcl2_first25])
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

    assert os.path.exists(f"{bcl2_first25}.hp.k24.scaled5.sig.zip.siglist")
    with open(f"{bcl2_first25}.hp.k24.scaled5.sig.zip.siglist") as f:
        assert f.readlines() == [f"{bcl2_first25}.hp.k24.scaled5.sig.zip"]

    # TODO: Actually test this, e.g. with what Sourmash Branchwater does for testing indices
    # https://github.com/sourmash-bio/sourmash_plugin_branchwater/blob/776727091d78af7de4990578ec6d07b12882c1b8/src/python/tests/test_index.py#L595
    # https://github.com/seanome/kmerseek/issues/3
    assert os.path.exists(f"{bcl2_first25}.hp.k24.scaled5.sig.zip.rocksdb")
