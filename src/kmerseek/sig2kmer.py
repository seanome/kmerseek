# ! /usr/bin/env python3
"""
Given a signature file and a collection of sequences, output all of the
k-mers and sequences that match a hashval in the signature file.

Cribbed from https://github.com/dib-lab/sourmash/pull/724/
"""

from typing import Literal
from tempfile import NamedTemporaryFile

import sourmash
import polars as pl
import screed


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

    def __init__(self, sig, fasta, moltype, ksize, scaled, save_kmers, save_sequences):
        self.signatures = [sig]  # List of signature files
        self.sequences = [fasta]  # List of sequence files
        self.ksize = ksize
        self.quiet = False
        self.force = True
        self.moltype = moltype
        self.translate = False
        self.check_sequence = True
        self.save_kmers = save_kmers
        self.save_sequences = save_sequences
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


def add_encoding_to_kmers_pl(kmers, moltype):

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

    return kmers_with_encoding


def add_start_position_to_kmers_pl(kmers, fasta):
    sequences = {}
    with screed.open(fasta) as records:
        for record in records:
            sequences[record["name"]] = record["sequence"]

    def find_start(row):
        sequence = sequences[row["sequence_name"]]
        kmer = row["kmer"]
        return sequence.index(kmer)

    # Add start column using map expression
    kmers_with_start = kmers.with_columns(
        pl.struct(["sequence_name", "kmer"])
        # Return unsigned integer32, range of 0 to 4,294,967,295
        # Don't expect protein sequences to be larger than 1 million really
        # position should always be positive, error if not
        .map_elements(find_start, return_dtype=pl.UInt32).alias("start")
    )
    return kmers_with_start


def postprocess_kmers(
    in_csv: str,
    in_fasta: str,
    moltype: Literal["hp", "dayhoff", "protein"],
    out_pq: str,
):
    """Add moltype-encoded version and start positions for each kmer

    Args:
        in_csv (str): _description_
        in_fasta (str): _description_
        moltype (Literal[&quot;hp&quot;, &quot;dayhoff&quot;, &quot;protein&quot;]): _description_
        out_pq (str): _description_
    """
    print(f"Reading in k-mers, adding {moltype} encoded values")
    # Read in as a LazyFrame to not load everything into memory
    kmers = pl.scan_csv(in_csv)

    kmers_with_encoding = add_encoding_to_kmers_pl(kmers, moltype)
    kmers_with_encoding_and_starts = add_start_position_to_kmers_pl(
        kmers_with_encoding, in_fasta
    )

    kmers_with_encoding_and_starts.sink_parquet(out_pq)


def get_kmers_cli(sig, fasta, moltype, ksize, scaled):
    # Create a temporary file for kmers CSV
    # TODO: Could potentially run out of temporary disk space but ...
    # maybe can stream the output to a parquet file?
    with NamedTemporaryFile(suffix=".csv") as tmp_kmers, NamedTemporaryFile(
        suffix=".fasta"
    ) as tmp_fasta:
        args = Args(
            sig=sig,
            fasta=fasta,
            moltype=moltype,  # or "dna", "dayhoff", "hp"
            ksize=ksize,
            scaled=scaled,
            save_kmers=tmp_kmers.name,  # Use the temp csv path for save_kmers
            save_sequences=tmp_fasta.name,  # Use the temp fasta path for save_sequences
        )

        # TODO: this "works" but calls what's normally a CLI as a Python method...
        # would love for this to be just a regular Python method that you could
        # feed sigs and fastas to
        print(f"Calling get_kmers_cli on {sig} with {fasta}")
        print(f"Saving matches to {args.save_kmers}")
        # CLI call: https://github.com/sourmash-bio/sourmash/blob/c209e7d39d80aa8eceed8d5a0c91568f96fded3f/src/sourmash/cli/sig/kmers.py#L96
        # Actual code that gets run: https://github.com/sourmash-bio/sourmash/blob/c209e7d39d80aa8eceed8d5a0c91568f96fded3f/src/sourmash/sig/__main__.py#L1087
        sourmash.sig.__main__.kmers(args)

        out_pq = _make_kmer_filename(sig)

        # Process the temp CSV file
        postprocess_kmers(tmp_kmers.name, tmp_fasta.name, moltype, out_pq)

        return out_pq
