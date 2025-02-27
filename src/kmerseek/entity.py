import os

import polars as pl

from .sketch import sketch, _make_sigfile
from .sig2kmer import get_kmers_cli, _make_kmer_filename
from .logging import logger


class KmerseekEntity:
    """Base class to be inherited by KmerseekQuery and KmerseekIndex"""

    def __init__(
        self, fasta, moltype, ksize, scaled, sketch=None, kmers_pq=None, force=False
    ):
        # These are all filenames of where the data is stored
        self.fasta = fasta
        self.sketch_kws = dict(moltype=moltype, ksize=ksize, scaled=scaled)
        self.force = force

        # If the sketch aka signature file is precomputed, then assign here
        if sketch is not None:
            self._sketch = sketch

        # If the kmers parquet file is precomputed, then assign here
        if kmers_pq is not None:
            self._kmers_pq = kmers_pq

    @property
    def sketch(self):
        """String of k-mer signature filename"""
        if not hasattr(self, "_sig") and not hasattr(self, "_precomputed_sig"):
            sigfile = _make_sigfile(self.fasta, **self.sketch_kws)
            if self.force or not os.path.exists(sigfile):
                if os.path.exists(sigfile):
                    logger.info(f"Found {sigfile} file, but re-making with '--force'")
                self._sketch = sketch(self.fasta, **self.sketch_kws)
            else:
                logger.info(
                    f"Found signature file {sigfile}, skipping! Re-make with '--force'"
                )
                self._sketch = sigfile
        return self._sketch

    @property
    def kmers_pq(self):
        """String of k-mer csv filename"""
        if not hasattr(self, "_kmers_pq"):
            pq = _make_kmer_filename(self.sketch)
            if self.force or not os.path.exists(pq):
                if os.path.exists(pq):
                    logger.info(f"Found {pq} file, but re-making with '--force'")
                self._kmers_pq = get_kmers_cli(
                    self.sketch, self.fasta, **self.sketch_kws
                )
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
