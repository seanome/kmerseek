import os

from click.testing import CliRunner


from kmerseek.main import cli


def test_search(ced9, bcl2_rocksdb):
    runner = CliRunner()
    result = runner.invoke(cli, ["search", "--ksize", "15", ced9, bcl2_rocksdb])
    assert result.exit_code == 0
