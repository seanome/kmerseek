import os

import polars as pl

from .logging import logger
from .sig2kmer import get_kmers_cli, _make_kmer_filename
from .sketch import sketch, _make_sigfile


class KmerseekEntity:
    """Base class to be inherited by KmerseekQuery and KmerseekIndex"""

    def __init__(self, fasta, moltype, ksize, scaled, force=False):
        # These are all filenames of where the data is stored
        self.fasta = fasta
        self.sketch_kws = dict(moltype=moltype, ksize=ksize, scaled=scaled)
        self.force = force

    @property
    def sig(self):
        """String of k-mer signature filename"""
        if not hasattr(self, "_sig"):
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
        """String of k-mer csv filename"""
        if not hasattr(self, "_kmers_pq"):
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
        """Polars LazyFrame of all the kmers"""
        if not hasattr(self, "_kmers"):
            self._kmers = pl.scan_parquet(self.kmers_pq)
        return self._kmers
