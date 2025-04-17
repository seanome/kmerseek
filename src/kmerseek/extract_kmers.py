import os
import polars as pl

from .logging import logger
from .sig2kmer import get_kmers_cli, _make_kmer_filename
from .sketch import sketch, _make_sigfile


class KmerseekExtractKmers:
    """Class dedicated to extracting k-mers from a FASTA file using sourmash."""

    def __init__(self, fasta, moltype, ksize, scaled, force=False):
        self.fasta = fasta
        self.sketch_kws = dict(moltype=moltype, ksize=ksize, scaled=scaled)
        self.force = force
        self._sig = None
        self._kmers_pq = None
        self._kmers_lazyframe = None

    @property
    def sig(self):
        """String of k-mer signature filename."""
        if self._sig is None:
            sigfile = _make_sigfile(self.fasta, **self.sketch_kws)
            if self.force or not os.path.exists(sigfile):
                if os.path.exists(sigfile):
                    logger.info(f"Found {sigfile} file, but re-making with '--force'")
                self._sig = sketch(self.fasta, **self.sketch_kws)
            else:
                logger.info(
                    f"Found signature file {sigfile}, skipping! Re-make with '--force'"
                )
                self._sig = sigfile
        return self._sig

    @property
    def kmers_pq(self):
        """String of k-mer parquet filename."""
        if self._kmers_pq is None:
            pq = _make_kmer_filename(self.sig)
            if self.force or not os.path.exists(pq):
                if os.path.exists(pq):
                    logger.info(f"Found {pq} file, but re-making with '--force'")
                self._kmers_pq = get_kmers_cli(self.sig, self.fasta, **self.sketch_kws)
            else:
                logger.info(
                    f"Found k-mer parquet {pq}, skipping! Re-make with '--force'"
                )
                self._kmers_pq = pq
        return self._kmers_pq

    @property
    def kmers_lazyframe(self):
        """Polars LazyFrame of all the kmers."""
        if self._kmers_lazyframe is None:
            if self.kmers_pq:
                self._kmers_lazyframe = pl.scan_parquet(self.kmers_pq)
            else:
                logger.warning("K-mer parquet file not generated.")
                return None
        return self._kmers_lazyframe

    def extract_kmers(self):
        """Explicitly trigger k-mer extraction and return the LazyFrame."""
        _ = self.kmers_pq  # Trigger k-mer generation if needed
        return self.kmers_lazyframe
