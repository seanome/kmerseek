from collections.abc import Iterable

import click
from sourmash_plugin_branchwater import sourmash_plugin_branchwater
import pandas as pd
import polars as pl

from .sig2kmer import _make_kmer_filename
from .uniprot import get_domains
from .sketch import make_sketch_kws
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


# TODO: benchmark manysearch vs multisearch: https://github.com/seanome/kmerseek/issues/2
def do_manysearch(query, target, output, ksize, scaled, moltype):
    query_siglist = _make_siglist_file(query.sig)

    # TODO: Figure out why target.rocksdb isn't working
    # https://github.com/seanome/kmerseek/issues/4
    target_siglist = _make_siglist_file(target.sig)
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


def do_multisearch(query, target, output, moltype, ksize, scaled):

    sourmash_plugin_branchwater.do_multisearch(
        query.sig,
        # TODO: Figure out why target.rocksdb isn't working
        # https://github.com/seanome/kmerseek/issues/4
        target.sig,
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
    merge_results_search_on = ["match_name", "query_name", "query_filename"]
    merge_results_kmers_on = [
        "sequence_name_match",
        "sequence_name_query",
        "sequence_file_query",
    ]

    merge_query_target_kmers_on = ["encoded", "hashval"]

    def __init__(
        self,
        output_csv: str,
        query: KmerseekQuery,
        target: KmerseekIndex,
    ):
        self.output_csv = output_csv
        self.query = query
        self.target = target
        self.results = self.read_results()

    def filter_target_kmers_in_results(self):
        # This might not actually be necessary since merge_query_target_kmers
        # would do the filter implicitly by only taking the inner matches
        self.target_kmers_in_results = self.target.kmers_lazyframe.filter(
            pl.col("sequence_name").is_in(self.results["match_name"])
        )

    def _prep_kmers_for_merging(self, kmers, suffix):
        cols_to_rename = ["kmer", "start", "sequence_name", "sequence_file"]

        renamer = {x: f"{x}{suffix}" for x in cols_to_rename}
        kmers_renamed = kmers.rename(renamer).sort(self.merge_query_target_kmers_on)
        return kmers_renamed

    def merge_query_target_kmers(self):

        query_kmers_prepped = self._prep_kmers_for_merging(
            self.query.kmers_lazyframe, "_query"
        )
        target_kmers_prepped = self._prep_kmers_for_merging(
            self.target.kmers_lazyframe, "_match"
        )

        import pdb

        pdb.set_trace()
        # No collect yet
        self.kmers = query_kmers_prepped.merge_sorted(
            target_kmers_prepped,
            key=self.merge_query_target_kmers_on,
        )

    def compute_tf_idf(self):
        raise NotImplementedError()

    def merge_results_kmers(self):

        self.results_with_kmers = self.results.join(
            self.kmers,
            left_on=self.merge_fastgather_on,
            right_on=self.merge_kmers_on,
        )

    def read_results(self):
        return pl.scan_csv(self.output_csv)

    def show_results_per_gene(self):
        results_per_gene = self.results_with_kmers.groupby("gene")

        show_results(results_per_gene)


# TODO: benchmark manysearch vs multisearch: https://github.com/seanome/kmerseek/issues/2
@click.command()
@click.argument("query_fasta")
@click.argument("target_fasta")
@click.option("--moltype", default="hp")
@click.option("--ksize", default=24)
@click.option("--scaled", default=5)
@click.option("--output", default="kmerseek_search.csv")
def search(
    query_fasta,
    target_fasta,
    moltype="hp",
    ksize=24,
    scaled=5,
    output="kmerseek_search.csv",
):

    sketch_kwargs = make_sketch_kws(moltype, ksize, scaled)

    query = KmerseekQuery(query_fasta, **sketch_kwargs)
    # Get kmers for the query
    _ = query.kmers_pq

    target = KmerseekIndex(target_fasta, **sketch_kwargs)

    do_multisearch(query, target, output, **sketch_kwargs)
    results = KmerseekResults(output, query, target)

    results.merge_query_target_kmers()
    results.merge_results_kmers()
