import csv
import os

import click
import sourmash_plugin_branchwater
import sourmash
import pandas as pd

from .sig2kmer import get_kmers
from .uniprot import get_domains
from .sketch import sketch


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
    stitched = ''
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



@click.command()
def main(
    query_fasta,
    target_sig,
    target_kmers,
    ksize=24,
    moltype="hp",
    scaled=5,
    search_output="search.csv",
):
    sketch_kwargs = dict(ksize=ksize, moltype=moltype, scaled=scaled)

    # Set defaults for search
    search_params = (
        dict(sketch_kwargs)
        .update("threshold", 0)
        .update("prob_significant_overlap", True)
        .update("output_path", search_output)
    )
    # Can this be async?
    query_sig = sketch(query_fasta, **sketch_kwargs)
    query_kmers = get_kmers(query_sig, query_fasta, **sketch_kwargs)

    sourmash_plugin_branchwater.do_multisearch(query_sig, target_sig, **search_params)

    results = pd.read_csv(search_output)

    results_with_kmers = results.join(query_kmers).join(target_kmers)

    results_per_gene = results_with_kmers.groupby("gene")

    show_results(results_per_gene)