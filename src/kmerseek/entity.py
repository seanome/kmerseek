import os

import polars as pl

from .sketch import sketch, _make_sigfile
from .sig2kmer import get_kmers_cli, _make_kmer_filename


class KmerseekEntity:
    """Base class to be inherited by KmerseekQuery and KmerseekIndex"""

    def __init__(self, fasta, moltype, ksize, scaled):
        # These are all filenames of where the data is stored
        self.fasta = fasta
        self.sketch_kws = dict(moltype=moltype, ksize=ksize, scaled=scaled)

    @property
    def sig(self):
        """String of k-mer signature filename"""
        if not hasattr(self, "_sig"):
            self._sig = sketch(self.fasta, **self.sketch_kws)
        return self._sig

    @property
    def kmers_csv(self):
        """String of k-mer csv filename"""
        if not hasattr(self, "_kmers_csv"):
            # kmers_csv = os.path.exists(_make_kmer_filename(self.sig))
            # if not kmers_csv:
            self._kmers_csv = get_kmers_cli(self.sig, self.fasta, **self.sketch_kws)
            # else:
            #     self._kmers_csv = kmers_csv
        return self._kmers_csv

    @property
    def kmers_lazyframe(self):
        """Polars LazyFrame of all the kmers"""
        # if not hasattr(self, "_kmers_csv"):
        #     self._kmers_csv = get_kmers(self.sig, self.fasta, **self.sketch_kws)
        if not hasattr(self, "_kmers"):
            self._kmers = pl.scan_csv(self._kmers_csv)
        return self._kmers
