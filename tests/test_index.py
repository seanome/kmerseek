import os

from click.testing import CliRunner


from kmerseek.main import cli


def test_index(bcl2_first25):
    runner = CliRunner()
    result = runner.invoke(cli, ["index", bcl2_first25])
    assert result.exit_code == 0

    # Make sure all the files got created
    assert os.path.exists(f"{bcl2_first25}.manysketch.csv")
    assert os.path.exists(f"{bcl2_first25}.hp.k24.scaled5.sig.zip")
    assert os.path.exists(f"{bcl2_first25}.hp.k24.scaled5.sig.zip.kmers.csv")
    assert os.path.exists(f"{bcl2_first25}.hp.k24.scaled5.sig.zip.siglist")
    assert os.path.exists(f"{bcl2_first25}.hp.k24.scaled5.sig.zip.rocksdb")
