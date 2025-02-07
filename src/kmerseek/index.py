from dataclasses import dataclass

import click
from sourmash.logging import notify
import pandas as pd

from sourmash_plugin_branchwater import sourmash_plugin_branchwater

from .entity import KmerseekEntity
from .sketch import sketch, make_sketch_kws
from .sig2kmer import get_kmers


class KmerseekIndex(KmerseekEntity):
    """Class to hold index locations for searching"""

    @property
    def rocksdb(self):
        if not hasattr(self, "_kmers_csv"):
            self._rocksdb = make_rocksdb_index(self.sig, **self.sketch_keywords)
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
    sketch_keywords = make_sketch_kws(moltype, ksize, scaled)

    kmerseek_index = KmerseekIndex(fasta, **sketch_keywords)
    sig = sketch(fasta, **sketch_keywords)

    kmers = get_kmers(sig, fasta, **sketch_keywords)

    rocksdb = make_rocksdb_index(sig, **sketch_keywords)
