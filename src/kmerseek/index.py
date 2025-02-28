import click
from sourmash.logging import notify
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
            self._rocksdb = make_rocksdb_index(self.sketch, **self.sketch_kws)
            # else:
            #     self._rocksdb = rocksdb
        return self._rocksdb


def _make_siglist_file(sig: str):
    siglist = f"{sig}.siglist"
    with open(siglist, "w") as f:
        f.write(f"{sig}")
    return siglist


def _make_rocksdb_filename(sig: str):
    return f"{sig}.rocksdb"


def make_rocksdb_index(sig, moltype, ksize, scaled):
    output = _make_rocksdb_filename(sig)
    notify(f"indexing all sketches in '{sig}'")

    siglist = _make_siglist_file(sig)

    # Colors set to false because that's what the sourmash_plugin_branchwater code does
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
        colors,
        internal_storage,
    )
    if status == 0:
        notify(f"...index is done! results in '{output}'")


@click.command()
@click.argument("fasta")
@click.option("--moltype", default="hp")
@click.option("--ksize", type=int, default=24)
@click.option("--scaled", type=int, default=5)
@click.option(
    "--force",
    is_flag=True,
    help="Force creation of signature, kmer parquet, and rocksdb even if they're already there",
)
@click.pass_context
def index(ctx, fasta, moltype="hp", ksize=24, scaled=5, force=False):
    """Prepare a database index for searching with k-mers"""
    # Now call each individual step
    ctx.forward(index_create_sketch)
    ctx.forward(index_create_kmers_pq)
    ctx.forward(index_create_rocksdb)


@click.command()
@click.argument("fasta")
@click.option("--moltype", default="hp")
@click.option("--ksize", type=int, default=24)
@click.option("--scaled", type=int, default=5)
@click.option(
    "--force",
    is_flag=True,
    help="Force creation of k-mer signature file even if file already exists",
)
def index_create_sketch(fasta, moltype="hp", ksize=24, scaled=5, force=False):
    """Substep of index: low memory, parallelized"""
    sketch_keywords = make_sketch_kws(moltype, ksize, scaled)

    kmerseek_index = KmerseekIndex(fasta, force=force, **sketch_keywords)
    _ = kmerseek_index.sketch


@click.command()
@click.argument("fasta")
@click.option("--moltype", default="hp")
@click.option("--ksize", type=int, default=24)
@click.option("--scaled", type=int, default=5)
@click.option(
    "--force",
    is_flag=True,
    help="Force creation of k-mer parquet file even if file already exists",
)
def index_create_kmers_pq(fasta, moltype="hp", ksize=24, scaled=5, force=False):
    """Substep of index: Extract k-mer sequences and encodings to a parquet file

    Low memory, may take a long time"""
    sketch_keywords = make_sketch_kws(moltype, ksize, scaled)

    kmerseek_index = KmerseekIndex(fasta, force=force, **sketch_keywords)
    _ = kmerseek_index.kmers_pq


@click.command()
@click.argument("fasta")
@click.option("--moltype", default="hp")
@click.option("--ksize", type=int, default=24)
@click.option("--scaled", type=int, default=5)
@click.option(
    "--force",
    is_flag=True,
    help="Force creation of rocksdb index even if already exists",
)
def index_create_rocksdb(fasta, moltype="hp", ksize=24, scaled=5, force=False):
    """Substep of index: Creates RocksDB index for fast searching"""
    sketch_keywords = make_sketch_kws(moltype, ksize, scaled)

    kmerseek_index = KmerseekIndex(fasta, force=force, **sketch_keywords)
    _ = kmerseek_index.rocksdb
