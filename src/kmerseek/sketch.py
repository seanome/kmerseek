
import sourmash_plugin_branchwater


def sketch(fasta, moltype, ksize, scaled):
    param_string = f"{moltype},k={ksize},scaled={scaled}"
    sigfile = f"{fasta}.{moltype}.k{ksize}.scaled{scaled}.sig.zip"
    sourmash_plugin_branchwater.do_manysketch(
        [fasta],
        param_string,
        output=sigfile,
        singleton=True,
        force=True,
    )
    return sigfile


