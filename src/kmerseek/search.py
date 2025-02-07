from collections.abc import Iterable

import click
from sourmash_plugin_branchwater import sourmash_plugin_branchwater
import pandas as pd

from .sig2kmer import get_kmers, _make_kmer_filename
from .uniprot import get_domains
from .sketch import sketch, make_sketch_kws
from .index import KmerseekIndex, _make_siglist_file
from .query import KmerseekQuery

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


def single_stitch_together_kmers(kmers: Iterable[str], i_kmers: Iterable[int]):
    stitched = ""
    for i, (i_kmer, kmer) in enumerate(zip(i_kmers, kmers)):
        if i == 0:
            stitched = kmer
            prev_i_kmer = i_kmer
        else:
            kmer_slice = i_kmer - prev_i_kmer
            stitched += kmer[-kmer_slice:]
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


def do_manysearch(query, target, output, ksize, scaled, moltype):
    query_siglist = _make_siglist_file(query.sig)
    target_siglist = _make_siglist_file(target.rocksdb)
    sourmash_plugin_branchwater.do_manysearch(
        query_siglist,
        target_siglist,
        0,  # threshold=0 to show all matches, even with only 1 k-mer
        ksize,
        scaled,
        moltype,
        output,
        False,  # Don't ignore abundance
        False,  # Don't output ALL comparisons, only the onse with at least a match
    )


def do_multisearch(query_sig, target_rocksdb, output, moltype, ksize, scaled):

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


class KmerseekResults:
    merge_fastgather_on = ["match_name", "query_name", "query_filename"]
    merge_kmers_on = ["read_name_match", "read_name_query", "filename_query"]

    def __init__(
        self,
        output_csv,
        query,
        target,
    ):
        self.output_csv = output_csv
        self.query = query
        self.target = target

    def merge_kmers(self):
        merge_on = ["kmer_in_alphabet", "hashval"]

        self.kmers = self.query.kmers.merge(
            self.target.kmers,
            left_on=merge_on,
            right_on=merge_on,
            suffixes=("_query", "_match"),
        )

    def compute_tf_idf(self):
        pass

    def merge_results_kmers(self):

        self.results_with_kmers = self.results.join(
            self.query_target_kmers,
            left_on=self.merge_fastgather_on,
            right_on=self.merge_kmers_on,
        )

    def read_results(self):
        self.results = pd.read_csv(self.output)

    def show_results_per_gene(self):
        results_per_gene = self.results_with_kmers.groupby("gene")

        show_results(results_per_gene)


@click.command()
@click.argument("query_fasta")
@click.argument("target_fasta")
@click.option("--moltype", default="hp")
@click.option("--ksize", default=24)
@click.option("--scaled", default=5)
@click.option("--output", default="kmerseek_search.csv")
@click.option(
    "--search-type", default="multisearch", help="either 'multisearch' or 'manysearch'"
)
def search(
    query_fasta,
    target_fasta,
    moltype="hp",
    ksize=24,
    scaled=5,
    output="kmerseek_search.csv",
    search_type="multisearch",
):

    sketch_kwargs = make_sketch_kws(moltype, ksize, scaled)

    query = KmerseekQuery(query_fasta, **sketch_kwargs)

    target = KmerseekIndex(target_fasta, **sketch_kwargs)

    if search_type == "multisearch":
        do_multisearch(query, target, output, **sketch_kwargs)
    elif search_type == "manysearch":
        do_manysearch(query, target, output, **sketch_kwargs)
