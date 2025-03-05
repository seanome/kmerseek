import click

from .index import (
    index,
    index_01_create_sketch,
    index_02_create_kmers_pq,
    index_03_create_rocksdb,
)
from .search import (
    search,
    search_01_create_query_sketch,
    search_02_create_query_kmers_pq,
    search_03_do_search,
    search_04_show_results,
)


@click.group()
def cli():
    """Kmerseek performs efficient protein domain annotation search with reduced amino acid k-mers"""
    pass


cli.add_command(index)
cli.add_command(index_01_create_sketch)
cli.add_command(index_02_create_kmers_pq)
cli.add_command(index_03_create_rocksdb)
cli.add_command(search)
cli.add_command(search_01_create_query_sketch)
cli.add_command(search_02_create_query_kmers_pq)
cli.add_command(search_03_do_search)
cli.add_command(search_04_show_results)
