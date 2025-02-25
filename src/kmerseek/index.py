import os

import click
from sourmash.logging import notify
import pandas as pd

from sourmash_plugin_branchwater import sourmash_plugin_branchwater

from .entity import KmerseekEntity
from .sketch import make_sketch_kws


class KmerseekIndex(KmerseekEntity):
    """Class to hold index locations for searching"""

    @property
    def rocksdb(self):
        if not hasattr(self, "_rocksdb"):
            # rocksdb = _make_rocksdb_filename(self.sig)
            # if not os.path.exists(rocksdb):
            self._rocksdb = make_rocksdb_index(self.sig, **self.sketch_kws)
            # else:
            #     self._rocksdb = rocksdb
        return self._rocksdb


def _make_siglist_file(sig):
    siglist = f"{sig}.siglist"
    with open(siglist, "w") as f:
        f.write(f"{sig}")
    return siglist


def _make_rocksdb_filename(sig):
    return f"{sig}.rocksdb"


def make_rocksdb_index(sig, moltype, ksize, scaled):
    output = _make_rocksdb_filename(sig)
    notify(f"indexing all sketches in '{sig}'")

    siglist = _make_siglist_file(sig)

    # Colros set to false because that's what the sourmash_plugin_branchwater code does
    colors = False

    # internal storage is false: don't store the signatures in the index, since the signatures are
    # in the same filesystem, "next door" in the same folder
    internal_storage = False

    status = sourmash_plugin_branchwater.do_index(
        siglist,
        ksize,
        scaled,
        moltype,
        output,
        colors,  # colors - currently must be false?
        internal_storage,
    )
    if status == 0:
        notify(f"...index is done! results in '{output}'")


@click.command()
@click.argument("fasta")
@click.option("--moltype", default="hp")
@click.option("--ksize", type=int, default=24)
@click.option("--scaled", type=int, default=5)
def index(fasta, moltype="hp", ksize=24, scaled=5):
    index_sketch(fasta, moltype, ksize, scaled)
    index_kmers_pq(fasta, moltype, ksize, scaled)
    index_rocksdb(fasta, moltype, ksize, scaled)


@click.command()
@click.argument("fasta")
@click.option("--moltype", default="hp")
@click.option("--ksize", type=int, default=24)
@click.option("--scaled", type=int, default=5)
def index_sketch(fasta, moltype="hp", ksize=24, scaled=5):
    sketch_keywords = make_sketch_kws(moltype, ksize, scaled)

    kmerseek_index = KmerseekIndex(fasta, **sketch_keywords)
    _ = kmerseek_index.sig


@click.command()
@click.argument("fasta")
@click.option("--moltype", default="hp")
@click.option("--ksize", type=int, default=24)
@click.option("--scaled", type=int, default=5)
def index_kmers_pq(fasta, moltype="hp", ksize=24, scaled=5):
    sketch_keywords = make_sketch_kws(moltype, ksize, scaled)

    kmerseek_index = KmerseekIndex(fasta, **sketch_keywords)
    _ = kmerseek_index.kmers_pq


@click.command()
@click.argument("fasta")
@click.option("--moltype", default="hp")
@click.option("--ksize", type=int, default=24)
@click.option("--scaled", type=int, default=5)
def index_rocksdb(fasta, moltype="hp", ksize=24, scaled=5):
    sketch_keywords = make_sketch_kws(moltype, ksize, scaled)

    kmerseek_index = KmerseekIndex(fasta, **sketch_keywords)
    _ = kmerseek_index.rocksdb
