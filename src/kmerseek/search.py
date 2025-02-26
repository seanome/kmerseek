import os
import sys
from tempfile import NamedTemporaryFile

import click
import polars as pl
from sourmash_plugin_branchwater import sourmash_plugin_branchwater

from .index import KmerseekIndex, _make_siglist_file
from .logging import setup_logging, logger
from .query import KmerseekQuery
from .sketch import make_sketch_kws, MOLTYPES
from .uniprot import get_domains


def show_results(results_per_gene):
    for gene, match in results_per_gene:
        domains = get_domains(gene, match.target_start, match.target_end)

        logger.info(
            f"Found {match.query_name}:{match.query_start}-{match.query_end} "
            f"in {match.target_name}:{match.target_start}-{match.target_end}"
        )
        if domains:
            for domain in domains:
                logger.info(f"Found: {domain.name} in {domain.start}-{domain.end}")
        logger.info(match.query_seq)
        logger.info(match.hp_kmer)
        logger.info(match.target_seq)


def single_stitch_together_kmers(kmers: pl.Series, i_kmers: pl.Series):
    stitched = ""
    logger.debug(f"Input kmers: {kmers}")
    logger.debug(f"Input i_kmers: {i_kmers}")

    for i, (i_kmer, kmer) in enumerate(zip(i_kmers, kmers)):
        logger.debug(f"Processing kmer {i}: position={i_kmer}, kmer={kmer}")
        if i == 0:
            stitched = kmer
            prev_i_kmer = i_kmer
            logger.debug(f"First kmer: {stitched}")
        else:
            kmer_slice = i_kmer - prev_i_kmer
            logger.debug(f"Slice size: {kmer_slice}")
            stitched += kmer[-kmer_slice:]
            logger.debug(f"After adding slice: {stitched}")
        prev_i_kmer = i_kmer

    logger.debug(f"Final stitched sequence: {stitched}")
    return stitched


def stitch_kmers_in_query_match_pair(
    df: pl.LazyFrame,
    kmer_match="kmer_match",
    kmer_alphabet="encoded",
    kmer_query="kmer_query",
    start_match="start_match",
    start_query="start_query",
) -> pl.DataFrame:
    logger.debug("\nProcessing new group:")
    logger.debug(f"Input dataframe:\n{df}")

    df = df.sort(start_query, descending=False)
    logger.debug(f"Sorted dataframe:\n{df}")

    match_name = df["match_name"][0]
    query_name = df["query_name"][0]

    query = single_stitch_together_kmers(df[kmer_query], df[start_match])
    alphabet = single_stitch_together_kmers(df[kmer_alphabet], df[start_query])
    match = single_stitch_together_kmers(df[kmer_match], df[start_match])

    logger.debug(f"query: {query}")
    logger.debug(f"alpha: {alphabet}")
    logger.debug(f"match: {match}")

    # All these lengths should be the same -> if not, something went wrong
    assert len(query) == len(alphabet)
    assert len(alphabet) == len(match)

    length = len(query)

    # 0-based, open intervals (same as Python indexing)
    match_start = df[start_match].min()
    query_start = df[start_query].min()
    match_end = match_start + length
    query_end = query_start + length

    to_print = (
        f"\n---\nQuery Name: {query_name}"
        f"\nMatch Name: {match_name}"
        f"\nquery: {query} ({query_start}-{query_end})\n"
        f"alpha: {alphabet}\n"
        f"match: {match} ({match_start}-{match_end})"
    )

    stitched_output = pl.DataFrame(
        {
            "match_name": [match_name],
            "query_name": [query_name],  # Add this line
            "query_start": [query_start],
            "query_end": [query_end],
            "query": [query],
            "match_start": [match_start],
            "match_end": [match_end],
            "match": [match],
            "to_print": [to_print],
        }
    )
    return stitched_output


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
        False,  # Don't output ALL comparisons, only the ones with at least a match
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
        False,  # Don't output ALL comparisons, only the ones with at least a match
        output,
    )


class KmerseekResults:
    join_results_search_on = ["match_name", "query_name"]
    join_results_kmers_on = [
        "sequence_name_match",
        "sequence_name_query",
    ]

    join_query_target_kmers_on = ["encoded", "hashval"]

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

    def read_results(self):
        df = pl.scan_csv(self.output_csv)
        return df

    def filter_target_kmers_in_results(self):
        # This might not actually be necessary since join_query_target_kmers
        # would do the filter implicitly by only taking the inner matches
        self.target_kmers_in_results = self.target.kmers_lazyframe.filter(
            pl.col("sequence_name").is_in(self.results["match_name"])
        )

    def _prep_kmers_for_merging(self, kmers, suffix):
        cols_to_rename = ["kmer", "start", "sequence_name", "sequence_file"]

        renamer = {x: f"{x}{suffix}" for x in cols_to_rename}
        kmers_renamed = kmers.rename(renamer).sort(self.join_query_target_kmers_on)
        return kmers_renamed

    def join_query_target_kmers(self):
        query_kmers_prepped = self._prep_kmers_for_merging(
            self.query.kmers_lazyframe, "_query"
        )
        target_kmers_prepped = self._prep_kmers_for_merging(
            self.target.kmers_lazyframe, "_match"
        )

        # No collect yet
        self.kmers = query_kmers_prepped.join(
            target_kmers_prepped, on=self.join_query_target_kmers_on
        )

    def compute_tf_idf(self):
        raise NotImplementedError()

    def join_results_kmers(self):
        self.results_with_kmers = self.results.join(
            self.kmers,
            left_on=self.join_results_search_on,
            right_on=self.join_results_kmers_on,
        )

    def stitch_kmers_per_gene(self):
        # Define the schema for the output of stitch_kmers_in_query_match_pair
        output_schema = {
            "match_name": pl.Utf8,
            "query_name": pl.Utf8,  # Add schema for query_name
            "query_start": pl.UInt32,
            "query_end": pl.UInt32,
            "query": pl.Utf8,
            "match_start": pl.UInt32,
            "match_end": pl.UInt32,
            "match": pl.Utf8,
            "to_print": pl.Utf8,
        }

        self.per_gene_stitched_kmers = (
            self.results_with_kmers.group_by("match_name")
            .map_groups(stitch_kmers_in_query_match_pair, schema=output_schema)
            .sort(["query_start", "query_end"])
        ).collect()

    def show_results_per_gene(self):
        print(
            self.per_gene_stitched_kmers.select(pl.col("to_print")).write_csv(
                quote_style="never", include_header=False
            ),
            file=sys.stderr,
        )


@click.command()
@click.argument("query_fasta")
@click.argument("target_fasta")
@click.option("--moltype", default="hp")
@click.option("--ksize", default=24)
@click.option("--scaled", default=5)
@click.option(
    "--output", default=None, help="If not specified, then output results to stdout"
)
@click.option(
    "--sourmash-search-csv",
    default=None,
    help=(
        "Store sourmash search results in this CSV. If not specified, then a temporary file is created. "
        "Mostly for debugging purposes"
    ),
)
@click.option("--debug", is_flag=True, help="Enable debug logging")
def search(
    query_fasta: str,
    target_fasta: str,
    moltype: MOLTYPES = "hp",
    ksize: int = 24,
    scaled: int = 5,
    output: str | None = None,
    sourmash_search_csv: str | None = None,
    debug: bool = False,
):
    """Search for k-mers in target sequences."""
    # Set up logging based on debug flag
    setup_logging(debug)

    logger.debug("Starting search with parameters:")
    logger.debug(f"Query FASTA: {query_fasta}")
    logger.debug(f"Target FASTA: {target_fasta}")
    logger.debug(f"Moltype: {moltype}")
    logger.debug(f"K-size: {ksize}")
    logger.debug(f"Scaled: {scaled}")
    logger.debug(f"Output: {output}")

    sketch_kwargs = make_sketch_kws(moltype, ksize, scaled)

    query = KmerseekQuery(query_fasta, **sketch_kwargs)
    _ = query.kmers_pq

    target = KmerseekIndex(target_fasta, **sketch_kwargs)

    temp_file = None

    csv_path, temp_file = create_sourmash_search_output(sourmash_search_csv, temp_file)

    try:
        # Run the search with the appropriate CSV file
        do_manysearch(query, target, csv_path, **sketch_kwargs)
        results = KmerseekResults(csv_path, query, target)

        # Process results
        results.join_query_target_kmers()
        results.join_results_kmers()
        results.stitch_kmers_per_gene()
        results.show_results_per_gene()

        results.per_gene_stitched_kmers.write_csv(output)

    finally:
        # Clean up temporary file if we created one
        if temp_file and os.path.exists(csv_path):
            os.unlink(csv_path)


def create_sourmash_search_output(sourmash_search_csv, temp_file):
    if sourmash_search_csv is None:
        # Create a temporary file that will be automatically cleaned up
        temp_file = NamedTemporaryFile(suffix=".csv", delete=False)
        csv_path = temp_file.name
        temp_file.close()  # Close but don't delete yet
    else:
        csv_path = sourmash_search_csv
    return csv_path, temp_file
