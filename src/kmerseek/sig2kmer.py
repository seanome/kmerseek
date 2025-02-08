# ! /usr/bin/env python3
"""
Given a signature file and a collection of sequences, output all of the
k-mers and sequences that match a hashval in the signature file.

Cribbed from https://github.com/dib-lab/sourmash/pull/724/
"""

from typing import Literal
import tempfile
from pathlib import Path

import sourmash
import polars as pl


def _make_kmer_filename(sig):
    return f"{sig}.kmers.pq"


# TODO: maybe add this to sourmash.sig.__main__.kmers because it's useful to us,
# could be useful to others
def encode_protein(sequence, moltype):
    """Convert protein sequence from 20-letter amino acid alphabet to degenerate alphabet"""
    # Pre-determine the alphabet function based on moltype
    if moltype == "hp":
        alphabet_func = sourmash._lowlevel.lib.sourmash_aa_to_hp
    elif moltype == "dayhoff":
        alphabet_func = sourmash._lowlevel.lib.sourmash_aa_to_dayhoff
    else:
        raise ValueError(f"Unknown moltype: {moltype}")

    # Convert the entire sequence to bytes once
    byte_sequence = sequence.encode("utf-8")

    # Apply the alphabet function to each byte in the sequence and join the result
    degenerate = b"".join(
        alphabet_func(letter.to_bytes(1, "big")) for letter in byte_sequence
    ).decode()

    return degenerate


class Args:
    """Dummy class to call what is normally a command line function from Python"""

    def __init__(self, sig, fasta, moltype, ksize, scaled, output):
        self.signatures = [sig]  # List of signature files
        self.sequences = [fasta]  # List of sequence files
        self.ksize = ksize
        self.quiet = False
        self.force = True
        self.moltype = moltype
        self.translate = False
        self.check_sequence = True
        self.save_kmers = output
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


def add_encoding_to_kmers_pl(
    in_csv: str, moltype: Literal["hp", "dayhoff", "protein"], out_pq: str
):
    print(f"Reading in k-mers, adding {moltype} encoded values")
    # Read in as a LazyFrame to not load everything into memory
    kmers = pl.scan_csv(in_csv)

    # Apply encoding using a UDF (User Defined Function)
    if moltype == "hp" or moltype == "dayhoff":
        kmers_with_encoding = kmers.with_columns(
            pl.col("kmer")
            # Need to specify return type as return_dtype=pl.Utf8, otherwise polars
            # refuses to sink to parquet because it doesn't know the datatype
            .map_elements(
                lambda x: encode_protein(x, moltype), return_dtype=pl.Utf8
            ).alias("encoded")
        )
    else:
        kmers_with_encoding = kmers

    kmers_with_encoding.sink_parquet(out_pq)


def get_kmers_cli(sig, fasta, moltype, ksize, scaled):
    # Create a temporary file for kmers CSV
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp_kmers:
        args = Args(
            sig=sig,
            fasta=fasta,
            moltype=moltype,  # or "dna", "dayhoff", "hp"
            ksize=ksize,
            scaled=scaled,
            output=tmp_kmers.name,  # Use the temp file path for save_kmers
        )

        # TODO: this "works" but calls what's normally a CLI as a Python method...
        # would love for this to be just a regular Python method that you could
        # feed sigs and fastas to
        print(f"Calling get_kmers_cli on {sig} with {fasta}")
        print(f"Saving matches to {args.save_kmers}")
        sourmash.sig.__main__.kmers(args)

        out_pq = _make_kmer_filename(sig)

        # Process the temp CSV file
        add_encoding_to_kmers_pl(tmp_kmers.name, moltype, out_pq)

        # Clean up the temporary file
        Path(tmp_kmers.name).unlink()
