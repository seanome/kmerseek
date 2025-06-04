import click

from .index import index
from .search import search


@click.group()
def cli():
    """Kmerseek performs efficient protein domain annotation search with reduced amino acid k-mers"""
    pass


cli.add_command(index)
cli.add_command(search)
