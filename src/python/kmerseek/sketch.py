import os
from typing import Literal

from sourmash_plugin_branchwater import sourmash_plugin_branchwater

MOLTYPES = Literal["hp", "dayhoff", "protein"]


def make_sketch_kws(moltype: MOLTYPES, ksize: int, scaled: int):
    return dict(ksize=ksize, moltype=moltype, scaled=scaled)


def _make_manysketch_csv(fasta: str):
    csv = f"{fasta}.manysketch.csv"
    basename = os.path.basename(fasta)
    with open(csv, "w") as f:
        f.write("name,genome_filename,protein_filename\n")
        # We're only indexing proteins sequences, so skip the genome (DNA) filenames
        f.write(f"{basename},,{fasta}\n")
    return csv


def _make_sigfile(fasta: str, moltype: MOLTYPES, ksize: int, scaled: int):
    sigfile = f"{fasta}.{moltype}.k{ksize}.scaled{scaled}.sig.zip"
    return sigfile


def sketch(fasta: str, moltype: MOLTYPES, ksize: int, scaled: int):
    param_string = f"{moltype},k={ksize},scaled={scaled},abund"
    sigfile = _make_sigfile(fasta, moltype, ksize, scaled)

    csv = _make_manysketch_csv(fasta)
    sourmash_plugin_branchwater.do_manysketch(
        csv,
        param_string,
        output=sigfile,
        singleton=True,
        force=True,
    )
    return sigfile
