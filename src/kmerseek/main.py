import click

from .index import (
    index,
    index_create_sketch,
    index_create_kmers_pq,
    index_create_rocksdb,
)
from .search import search


@click.group()
def cli():
    """Kmerseek performs efficient protein domain annotation search with reduced amino acid k-mers"""
    pass


cli.add_command(index)
cli.add_command(index_create_sketch)
cli.add_command(index_create_kmers_pq)
cli.add_command(index_create_rocksdb)
cli.add_command(search)
