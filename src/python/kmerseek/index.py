import click
from sourmash.logging import notify
from sourmash_plugin_branchwater import sourmash_plugin_branchwater

from .entity import KmerseekEntity
from .logging import setup_logging, logger
from .sketch import make_sketch_kws


class KmerseekIndexBase(KmerseekEntity):
    """Base class for Kmerseek index functionality."""

    @property
    def rocksdb(self):
        if not hasattr(self, "_rocksdb"):
            self._rocksdb = make_rocksdb_index(self.sig, **self.sketch_kws)
        return self._rocksdb

    def __init__(self, fasta, moltype, ksize, scaled, force=False):
        super().__init__(
            fasta, moltype, ksize, scaled, force=force, extract_kmers=False
        )


class KmerseekIndexWithKmerExtraction(KmerseekIndexBase):
    """Kmerseek index that also performs k-mer extraction."""

    def __init__(self, fasta, moltype, ksize, scaled, force=False):
        super().__init__(fasta, moltype, ksize, scaled, force=force)
        self.extract_kmers = True
        # Trigger k-mer extraction upon initialization
        _ = self.kmers_pq
        _ = self.kmers_lazyframe


class KmerseekIndexWithoutKmerExtraction(KmerseekIndexBase):
    """Kmerseek index that does not perform k-mer extraction."""

    def __init__(self, fasta, moltype, ksize, scaled, force=False):
        super().__init__(fasta, moltype, ksize, scaled, force=force)
        self.extract_kmers = False


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
@click.option("--extract-kmers", is_flag=True, default=False)
@click.option("--debug", is_flag=True, help="Enable debug logging")
@click.option(
    "--force",
    is_flag=True,
    help="Force creation of signature, kmer parquet, and rocksdb even if they're already there",
)
def index(
    fasta,
    moltype="hp",
    ksize=24,
    scaled=5,
    extract_kmers=False,
    debug=False,
    force=False,
):
    # Set up logging based on debug flag
    setup_logging(debug)

    sketch_keywords = make_sketch_kws(moltype, ksize, scaled)

    if extract_kmers:
        kmerseek_index = KmerseekIndexWithKmerExtraction(
            fasta, force=force, **sketch_keywords
        )
        logger.info("K-mer extraction will be performed during indexing.")
    else:
        kmerseek_index = KmerseekIndexWithoutKmerExtraction(
            fasta, force=force, **sketch_keywords
        )
        logger.info("K-mer extraction will be skipped during indexing.")

        # Trigger the creation of signature and rocksdb for both cases
    _ = kmerseek_index.sig
    _ = kmerseek_index.rocksdb
    if extract_kmers:
        logger.info(f"K-mers stored in: {kmerseek_index.kmers_pq}")
