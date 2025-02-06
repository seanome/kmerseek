import os

from sourmash_plugin_branchwater import sourmash_plugin_branchwater


def _make_manysketch_csv(fasta):
    csv = f"{fasta}.manysketch.csv"
    basename = os.path.basename(fasta)
    with open(csv, "w") as f:
        f.write("name,genome_filename,protein_filename\n")
        # We're only indexing proteins sequences, so skip the genome (DNA) filenames
        f.write(f"{basename},,{fasta}\n")
    return csv


def sketch(fasta, moltype, ksize, scaled):
    param_string = f"{moltype},k={ksize},scaled={scaled}"
    sigfile = f"{fasta}.{moltype}.k{ksize}.scaled{scaled}.sig.zip"

    csv = _make_manysketch_csv(fasta)
    sourmash_plugin_branchwater.do_manysketch(
        csv,
        param_string,
        output=sigfile,
        singleton=True,
        force=True,
    )
    return sigfile
