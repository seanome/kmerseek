from click.testing import CliRunner
import pytest


from kmerseek.index import index



def test_hello_world(bcl2_first25):
  runner = CliRunner()
  result = runner.invoke(index, [bcl2_first25])
  assert result.exit_code == 0
  assert result.output == 'Hello Peter!\n'