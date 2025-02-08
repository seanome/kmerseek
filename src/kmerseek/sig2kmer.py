# ! /usr/bin/env python3
"""
Given a signature file and a collection of sequences, output all of the
k-mers and sequences that match a hashval in the signature file.

Cribbed from https://github.com/dib-lab/sourmash/pull/724/
"""

import sourmash


class Args:
    """Dummy class to call what is normally a command line function from Python"""

    def __init__(self, sig, fasta, moltype, ksize, scaled):
        self.signatures = [sig]  # List of signature files
        self.sequences = [fasta]  # List of sequence files
        self.ksize = ksize
        self.quiet = False
        self.force = True
        self.moltype = moltype
        self.translate = False
        self.check_sequence = True
        self.save_kmers = None
        self.save_sequences = None
        self.picklist = None
        self.scaled = scaled
        self.dna = False
        self.skipm1n3 = False
        self.skipm2n3 = False
        self.from_file = False

        self.dayhoff = False
        self.hp = False
        self.protein = False
        if moltype == "protein":
            self.protein = True
        elif moltype == "dayhoff":
            self.dayhoff = True
        elif moltype == "hp":
            self.hp = True


def get_kmers_cli(sig, fasta, moltype, ksize, scaled):
    # Example usage:
    args = Args(
        sig=sig,
        fasta=fasta,
        moltype=moltype,  # or "dna", "dayhoff", "hp"
        ksize=ksize,
        scaled=scaled,
    )

    # TODO: this "works" but calls what's normally a CLI as a Python method...
    # would love for this to be just a regular Python method that you could
    # feed sigs and fastas to
    sourmash.sig.__main__.kmers(args)
