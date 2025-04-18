import os

from .extract_kmers import KmerseekKmerExtractor  # Import the new class
from .logging import logger
from .sketch import sketch, _make_sigfile


class KmerseekEntity:
    """Base class to be inherited by KmerseekQuery and KmerseekIndex"""

    def __init__(self, fasta, moltype, ksize, scaled, force=False, extract_kmers=False):
        # These are all filenames of where the data is stored
        self.fasta = fasta
        self.sketch_kws = dict(moltype=moltype, ksize=ksize, scaled=scaled)
        self.force = force
        self.extract_kmers = extract_kmers
        self._sig = None
        self._kmer_extractor = None  # Instance of KmerseekKmerExtractor

    @property
    def sig(self):
        """String of k-mer signature filename"""
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
    def kmer_extractor(self):
        """Instance of KmerseekKmerExtractor."""
        if self.extract_kmers and self._kmer_extractor is None:
            self._kmer_extractor = KmerseekKmerExtractor(
                self.fasta,
                self.sketch_kws["moltype"],
                self.sketch_kws["ksize"],
                self.sketch_kws["scaled"],
                self.force,
            )
        return self._kmer_extractor

    @property
    def kmers_pq(self):
        """String of k-mer parquet filename"""
        if self.extract_kmers:
            if self.kmer_extractor:
                return self.kmer_extractor.kmers_pq
            else:
                return None
        else:
            logger.info("Skipping k-mer extraction")
            return None

    @property
    def kmers_lazyframe(self):
        """Polars LazyFrame of all the kmers"""
        if self.extract_kmers:
            if self.kmer_extractor:
                return self.kmer_extractor.extract_kmers()
            else:
                return None
        else:
            logger.info("K-mer extraction was skipped during initialization.")
            return None
