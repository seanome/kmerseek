import click
from sourmash_plugin_branchwater import sourmash_plugin_branchwater
import pandas as pd

from .sig2kmer import get_kmers
from .uniprot import get_domains
from .sketch import sketch, make_sketch_kws


protein_moltypes = "protein", "dayhoff", "hp"


def show_results(results_per_gene):
    for gene, match in results_per_gene:
        # This function currently doesn't do anything, just here for
        domains = get_domains(gene, match.target_start, match.target_end)

        print(
            f"Found {match.query_name}:{match.query_start}-{match.query_end} "
            f"in {match.target_name}:{match.target_start}-{match.target_end}"
        )
        if domains:
            for domain in domains:
                print(f"Found: {domain.name} in {domain.start}-{domain.end}")
        print(match.query_seq)
        print(match.hp_kmer)
        print(match.target_seq)


def single_stitch_together_kmers(kmer_series):
    stitched = ""
    for i, kmer in enumerate(kmer_series):
        if i == 0:
            stitched = kmer
        else:
            stitched += kmer[-1]
    return stitched


def stitch_together_species_hp_kmers(
    df,
    kmer_query="kmer_query",
    kmer_hp="kmer_hp",
    kmer_found="kmer_found",
    i_query="i_query",
):
    df = df.sort_values(i_query)
    query_stitched = single_stitch_together_kmers(df[kmer_query])
    hp_stitched = single_stitch_together_kmers(df[kmer_hp])
    found_stitched = single_stitch_together_kmers(df[kmer_found])
    start = df.i_query.min()
    end = start + len(found_stitched)
    to_print = (
        f"botryllus: {query_stitched}\nhp: {hp_stitched}\nhuman: {found_stitched}"
    )
    return pd.Series([start, end, to_print], index=["start", "end", "matches"])


def search(query_sig, target_rocksdb, output, moltype, ksize, scaled):

    sourmash_plugin_branchwater.do_multisearch(
        query_sig,
        target_rocksdb,
        0,  # threshold=0 to show all matches, even with only 1 k-mer
        ksize,
        scaled,
        moltype,
        False,  # No Average nucleotide identity (ANI) calculation
        True,  # Yes, calculate probability of overlap between target and query
        False,  # Don't output ALL comparisons, only the onse with at least a match
        output,
    )


@click.command()
@click.argument("query_fasta")
@click.argument("target_rocksdb")
@click.argument("target_kmers")
@click.option("--moltype", default="hp")
@click.option("--ksize", default=24)
@click.option("--scaled", default=5)
@click.option("--output", default="search.csv")
def search(
    query_fasta,
    target_rocksdb,
    target_kmers,
    moltype="hp",
    ksize=24,
    scaled=5,
    output="search.csv",
):
    sketch_kwargs = make_sketch_kws(moltype, ksize, scaled)

    query_sig = sketch(query_fasta, **sketch_kwargs)

    # Can this be async relative to the search? It can happen simultaneously
    query_kmers = get_kmers(query_sig, query_fasta, **sketch_kwargs)

    search(query_sig, target_rocksdb, output, **sketch_kwargs)

    results = pd.read_csv(output)

    results_with_kmers = results.join(query_kmers).join(target_kmers)

    results_per_gene = results_with_kmers.groupby("gene")

    show_results(results_per_gene)
