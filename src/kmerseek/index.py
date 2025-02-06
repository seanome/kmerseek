

import click
from sourmash.logging import notify

import sourmash_plugin_branchwater

from .sketch import sketch
from .sig2kmer import get_kmers

def make_rocksdb_index(sig, moltype, ksize, scaled):
    output = f'{sig}.rocksdb'
    notify(f"indexing all sketches in '{sig}'")
    status = sourmash_plugin_branchwater.do_index(
        [sig],
        ksize,
        scaled,
        moltype,
        output,
        False,  # colors - currently must be false?
        internal_storage=False
    )
    if status == 0:
        notify(f"...index is done! results in '{output}'")

@click.command()
@click.argument('fasta', help='input sequence file')
@click.option('--moltype')
@click.option('--ksize', type=int)
@click.option('--scaled', type=int)
def index(fasta, moltype='hp', ksize=24, scaled=5):
    sketch_keywords = dict(moltype='hp', ksize=24, scale=5)
    sig = sketch(fasta, **sketch_keywords)

    kmers = get_kmers(sig, fasta, **sketch_keywords)

    rocksdb = make_rocksdb_index(sig)

