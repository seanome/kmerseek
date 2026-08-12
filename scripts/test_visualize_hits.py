"""Tests for visualize_hits.py.

Run with: /Users/olga/anaconda3/envs/kmerseek-dev/bin/python3 -m pytest scripts/test_visualize_hits.py -v
"""

import csv
import os
import pathlib
import re
import sys

import matplotlib
import pytest

sys.path.insert(0, os.path.dirname(__file__))
import visualize_hits as vh

TESTDATA_DIR = os.path.join(os.path.dirname(__file__), "..", "tests", "testdata")
CED9_FASTA = os.path.join(TESTDATA_DIR, "fasta", "ced9.fasta")
BCL2_FASTA = os.path.join(TESTDATA_DIR, "fasta", "bcl2.fasta")
MULTI_FASTA = os.path.join(
    TESTDATA_DIR, "index",
    "bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta",
)

CED9_NAME = "sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1"
CED9_LEN = 280
BCL2_NAME = "sp|P10415|BCL2_HUMAN Apoptosis regulator Bcl-2 OS=Homo sapiens OX=9606 GN=BCL2 PE=1 SV=2"
BCL2_LEN = 239
SEPT4_NAME = "sp|O43236|SEPT4_HUMAN Septin-4 OS=Homo sapiens OX=9606 GN=SEPTIN4 PE=1 SV=1"
SEPT4_LEN = 478
B2L11_NAME = "sp|O43521|B2L11_HUMAN Bcl-2-like protein 11 OS=Homo sapiens OX=9606 GN=BCL2L11 PE=1 SV=1"
B2L11_LEN = 198


# --- read_fasta_lengths ---------------------------------------------------

def test_read_fasta_lengths_single_record_fasta():
    assert vh.read_fasta_lengths(CED9_FASTA) == {CED9_NAME: CED9_LEN}


def test_read_fasta_lengths_multi_record_fasta():
    lengths = vh.read_fasta_lengths(MULTI_FASTA)
    assert len(lengths) == 25
    assert lengths[SEPT4_NAME] == SEPT4_LEN
    assert lengths[B2L11_NAME] == B2L11_LEN


def test_read_fasta_lengths_no_records():
    assert vh.read_fasta_lengths(os.devnull) == {}


# --- short_label -----------------------------------------------------------

def test_short_label_uniprot_header():
    assert vh.short_label(BCL2_NAME) == "BCL2_HUMAN"


def test_short_label_no_pipes_uses_first_word():
    assert vh.short_label("myprotein some description") == "myprotein"


def test_short_label_truncates_long_names():
    name = "sp|P00000|" + "X" * 40 + " description"
    result = vh.short_label(name, max_len=10)
    assert result == "XXXXXXXXX…"
    assert len(result) == 10


# --- safe_filename -----------------------------------------------------------

def test_safe_filename_strips_special_characters():
    assert vh.safe_filename(BCL2_NAME) == "BCL2_HUMAN"


def test_safe_filename_no_pipes_uses_first_word():
    assert vh.safe_filename("weird!!@@## name") == "weird"


# --- resolve_query_names -----------------------------------------------------

def test_resolve_query_names_none_returns_all_sorted():
    names = {SEPT4_NAME, B2L11_NAME}
    assert vh.resolve_query_names(None, names) == sorted(names)


def test_resolve_query_names_exact_header_match():
    names = {CED9_NAME, BCL2_NAME}
    assert vh.resolve_query_names(CED9_NAME, names) == [CED9_NAME]


def test_resolve_query_names_short_gene_symbol_matches():
    names = {CED9_NAME, BCL2_NAME}
    assert vh.resolve_query_names("CED9", names) == [CED9_NAME]


def test_resolve_query_names_is_case_insensitive():
    names = {CED9_NAME, BCL2_NAME}
    assert vh.resolve_query_names("ced9_caeel", names) == [CED9_NAME]


def test_resolve_query_names_no_match_returns_empty():
    names = {CED9_NAME, BCL2_NAME}
    assert vh.resolve_query_names("NOT_A_GENE", names) == []


# --- merge_regions_by_target / _build_hit -----------------------------------

# region_subseq/target_subseq below are real residues 1-10 of CED9 and BCL2
# (true Bcl-2-family homologs), and moltype_seq is their actual thomas_dill
# HP encoding (h: A,C,F,I,L,M,V,W,Y; p: D,E,G,H,K,N,P,Q,R,S,T) -- not a made-up motif.
#
# The keys here are a contract with SearchResultCsv in src/rust/search.rs: this
# fixture is what the tests see instead of real `kmerseek search` output, so if it
# drifts from the struct the whole suite can pass while the script is broken against
# every real CSV. test_fixture_matches_rust_csv_schema below pins the two together.
#
# Values are native int/float/str, not strings-of-numbers: scan_csv() infers a
# schema, so iter_query_rows()'s to_dicts() hands the rest of this module typed
# values, never the all-string rows a stdlib csv.DictReader would produce.
def _row(target_name="tgt", region_start=0, region_end=10, region_length=10,
         containment=0.5, jaccard=0.1, query_enrichment=1.0, query_poisson_pvalue=0.05,
         region_poisson_pvalue=0.01, region_enrichment=2.0,
         moltype="hp_thomas_dill", region_subseq="MTRCTADNSL",
         moltype_seq="hpphphppph", target_subseq="MAHAGRTGYD", query_name=CED9_NAME):
    return {
        "query_name": query_name,
        "query_md5": "qmd5",
        "target_name": target_name,
        "target_md5": "tmd5",
        "containment": containment,
        "n_intersecting_hashes": 5,
        "ksize": 10,
        "scaled": 1,
        "moltype": moltype,
        "jaccard": jaccard,
        "max_containment": containment,
        "average_abund": 1.0,
        "median_abund": 1.0,
        "std_abund": 0.0,
        "containment_target_in_query": containment,
        "f_weighted_target_in_query": 1.0,
        "query_tfidf": 1.0,
        "mean_matched_kmer_freq": 0.1,
        "sum_matched_kmer_freq": 0.5,
        "query_expected_shared_kmers": 1.0,
        "query_enrichment": query_enrichment,
        "joint_kmer_freq": 0.0,
        "query_poisson_pvalue": query_poisson_pvalue,
        "region_search_space": 271,
        "db_n_targets": 25,
        "db_n_kmers": 7629,
        "run_n_queries": 1,
        "region_start": region_start,
        "region_end": region_end,
        "region_subseq": region_subseq,
        "target_start": region_start,
        "target_end": region_end,
        "target_subseq": target_subseq,
        "moltype_seq": moltype_seq,
        "region_length": region_length,
        "region_n_shared_kmers": 1,
        "region_expected_shared_kmers": 0.68,
        "region_poisson_pvalue": region_poisson_pvalue,
        "region_enrichment": region_enrichment,
    }


def test_fixture_matches_rust_csv_schema():
    """_row() must carry exactly the columns `kmerseek search` writes.

    Renaming a column in SearchResultCsv without updating this fixture is
    otherwise invisible: every test here keeps passing against the stale shape
    while the script raises KeyError on real output.
    """
    search_rs = pathlib.Path(__file__).resolve().parent.parent / "src" / "rust" / "search.rs"
    source = search_rs.read_text()
    struct_body = re.search(r"pub struct SearchResultCsv \{(.*?)\n\}", source, re.S).group(1)
    rust_columns = set(re.findall(r"^\s*pub (\w+):", struct_body, re.M))

    assert rust_columns, "failed to parse SearchResultCsv fields"
    assert set(_row()) == rust_columns


def test_merge_regions_within_gap_becomes_one_hit():
    rows = [
        _row(region_start=0, region_end=10, region_length=10),
        _row(region_start=15, region_end=25, region_length=10),  # gap of 5
    ]
    hits = vh.merge_regions_by_target(rows, gap_merge=10)
    assert len(hits) == 1
    assert hits[0]["start"] == 0
    assert hits[0]["end"] == 25
    assert hits[0]["n_regions"] == 2


def test_merge_regions_beyond_gap_stays_separate():
    rows = [
        _row(region_start=0, region_end=10, region_length=10),
        _row(region_start=25, region_end=35, region_length=10),  # gap of 15
    ]
    hits = vh.merge_regions_by_target(rows, gap_merge=10)
    assert len(hits) == 2
    assert (hits[0]["start"], hits[0]["end"]) == (0, 10)
    assert (hits[1]["start"], hits[1]["end"]) == (25, 35)


def test_merge_regions_groups_by_target_independently():
    rows = [
        _row(target_name="A", region_start=0, region_end=10, region_length=10),
        _row(target_name="B", region_start=5, region_end=15, region_length=10),
    ]
    hits = vh.merge_regions_by_target(rows, gap_merge=10)
    assert {h["target_name"] for h in hits} == {"A", "B"}
    assert len(hits) == 2


def test_build_hit_region_rows_keeps_every_region_sorted_by_start():
    # Every region must stay visible (not just one "representative" pick) so a
    # multi-region hit like BAK_HUMAN's two separate matches both get shown.
    # CED9 residues 21-25 ("ATGEM") and 1-20 ("MTRCTADNSLTNPAYRRRTM"), aligned
    # against the corresponding BCL2 residues -- two real, disjoint matches.
    second_region = _row(region_start=20, region_end=25, region_length=5,
                          region_subseq="ATGEM", moltype_seq="hppph",
                          target_subseq="YKLSQ", containment=0.9)
    first_region = _row(region_start=0, region_end=20, region_length=20,
                         region_subseq="MTRCTADNSLTNPAYRRRTM",
                         moltype_seq="hpphphppphppphhpppph",
                         target_subseq="NREIVMKYIHYKLSQRGYEW", containment=0.3)
    hit = vh._build_hit("tgt", [second_region, first_region])
    assert [row["region_subseq"] for row in hit["region_rows"]] == [
        "MTRCTADNSLTNPAYRRRTM", "ATGEM",
    ]


def test_build_hit_containment_is_max_across_cluster():
    rows = [_row(containment=0.2), _row(containment=0.7), _row(containment=0.4)]
    hit = vh._build_hit("tgt", rows)
    assert hit["containment"] == 0.7


def test_build_hit_carries_jaccard_enrichment_pvalue():
    # These are query-target *result* stats (same value on every region row of
    # a hit, since it's the same query-target pair), not per-region stats.
    rows = [_row(jaccard=0.123, query_enrichment=4.5, query_poisson_pvalue=1e-06)]
    hit = vh._build_hit("tgt", rows)
    assert hit["jaccard"] == 0.123
    assert hit["query_enrichment"] == 4.5
    assert hit["query_poisson_pvalue"] == 1e-06


# --- benjamini_hochberg / _target_pvalues -------------------------------------

def test_target_pvalues_one_entry_per_distinct_target():
    rows = [_row(target_name="A", region_poisson_pvalue=0.01),
            _row(target_name="A", region_poisson_pvalue=0.01),
            _row(target_name="B", region_poisson_pvalue=0.02)]
    assert vh._target_pvalues(rows) == {"A": 0.01, "B": 0.02}


def test_target_pvalues_keeps_strongest_region_per_target():
    # Region p-values differ row to row (unlike the query-level stats), so a target
    # is represented by its best region, not an arbitrary one.
    rows = [_row(target_name="A", region_poisson_pvalue=0.4),
            _row(target_name="A", region_poisson_pvalue=1e-08),
            _row(target_name="A", region_poisson_pvalue=0.2)]
    assert vh._target_pvalues(rows) == {"A": 1e-08}


def test_target_pvalues_corrects_region_scope_not_query_scope():
    # The case region scoring exists for: a diluted whole-query p-value next to a
    # strong region. Correcting the query scope would bury it at q~1.
    rows = [_row(target_name="cryptic", query_poisson_pvalue=0.99, region_poisson_pvalue=0.0007)]
    assert vh._target_pvalues(rows) == {"cryptic": 0.0007}


def test_benjamini_hochberg_single_pvalue_is_unchanged():
    assert vh.benjamini_hochberg({"x": 0.5}) == {"x": 0.5}


def test_benjamini_hochberg_known_example():
    # p.adjust(c(d=0.005, a=0.01, c=0.03, b=0.04), method="BH") in R gives
    # d=0.02, a=0.02, c=0.04, b=0.04.
    pvalues = {"a": 0.01, "b": 0.04, "c": 0.03, "d": 0.005}
    q = vh.benjamini_hochberg(pvalues)
    assert q["d"] == pytest.approx(0.02)
    assert q["a"] == pytest.approx(0.02)
    assert q["c"] == pytest.approx(0.04)
    assert q["b"] == pytest.approx(0.04)


def test_benjamini_hochberg_never_exceeds_one():
    q = vh.benjamini_hochberg({"x": 0.9, "y": 0.95})
    assert all(v <= 1.0 for v in q.values())


# --- _union_coverage ------------------------------------------------------

def test_union_coverage_non_overlapping_regions_sums_lengths():
    assert vh._union_coverage([(0, 10), (20, 30)]) == 20


def test_union_coverage_overlapping_regions_counts_union_once():
    # Heavily overlapping sliding-window matches must not be double-counted:
    # (0, 10) and (5, 15) share residues 5-10, so the union is 0-15 = 15, not 20.
    assert vh._union_coverage([(0, 10), (5, 15)]) == 15


def test_union_coverage_identical_regions_counts_once():
    assert vh._union_coverage([(0, 10), (0, 10), (0, 10)]) == 10


def test_build_hit_coverage_is_union_not_sum():
    rows = [
        _row(region_start=0, region_end=10, region_length=10),
        _row(region_start=5, region_end=15, region_length=10),
    ]
    hit = vh._build_hit("tgt", rows)
    assert hit["coverage"] == 15
    assert hit["end"] - hit["start"] == 15


# --- get_palette ------------------------------------------------------------

def test_get_palette_small_n_uses_set2():
    palette = vh.get_palette(4)
    assert palette == list(matplotlib.colormaps["Set2"].colors)


def test_get_palette_ten_targets_uses_tab10():
    palette = vh.get_palette(10)
    assert palette == list(matplotlib.colormaps["tab10"].colors)


def test_get_palette_twelve_targets_uses_set3():
    palette = vh.get_palette(12)
    assert palette == list(matplotlib.colormaps["Set3"].colors)


def test_get_palette_twenty_targets_uses_tab20():
    palette = vh.get_palette(20)
    assert palette == list(matplotlib.colormaps["tab20"].colors)


def test_get_palette_more_than_twenty_falls_back_to_tab20():
    palette = vh.get_palette(35)
    assert palette == list(matplotlib.colormaps["tab20"].colors)


# --- assign_lanes ------------------------------------------------------------

def test_assign_lanes_non_overlapping_hits_share_lane():
    hits = [{"start": 0, "end": 10}, {"start": 20, "end": 30}]
    n_lanes = vh.assign_lanes(hits)
    assert n_lanes == 1
    assert hits[0]["lane"] == 0
    assert hits[1]["lane"] == 0


def test_assign_lanes_overlapping_hits_get_separate_lanes():
    hits = [{"start": 0, "end": 20}, {"start": 10, "end": 30}]
    n_lanes = vh.assign_lanes(hits)
    assert n_lanes == 2
    lanes = {(h["start"], h["end"]): h["lane"] for h in hits}
    assert lanes[(0, 20)] != lanes[(10, 30)]


def test_assign_lanes_three_way_overlap_needs_three_lanes():
    hits = [{"start": 0, "end": 30}, {"start": 5, "end": 25}, {"start": 10, "end": 20}]
    assert vh.assign_lanes(hits) == 3


# --- wrap_seq ------------------------------------------------------------

def test_wrap_seq_short_sequence_stays_one_line():
    lines = vh.wrap_seq("MKVLLL")
    assert lines == ["MKVLLL"]


def test_wrap_seq_wraps_at_configured_width():
    seq = "A" * (vh.SEQ_WRAP_WIDTH + 5)
    lines = vh.wrap_seq(seq)
    assert len(lines) == 2
    assert lines[0] == "A" * vh.SEQ_WRAP_WIDTH
    assert lines[1] == "A" * 5


def test_wrap_seq_truncates_past_max_lines():
    seq = "A" * (vh.SEQ_WRAP_WIDTH * (vh.MAX_SEQ_LINES + 2))
    lines = vh.wrap_seq(seq)
    assert len(lines) == vh.MAX_SEQ_LINES
    assert lines[-1].endswith("…")


# --- plot_gene (integration smoke test) --------------------------------------

def test_plot_gene_writes_png_and_svg(tmp_path):
    hits = vh.merge_regions_by_target(
        [_row(target_name="TGT_A", region_start=10, region_end=40, region_length=30)],
        gap_merge=10,
    )
    png_path = tmp_path / "gene.hits.png"
    svg_path = tmp_path / "gene.hits.svg"
    vh.plot_gene(CED9_NAME, CED9_LEN, hits, [str(png_path), str(svg_path)])
    assert png_path.exists() and png_path.stat().st_size > 0
    assert svg_path.exists() and svg_path.stat().st_size > 0


def test_plot_gene_svg_text_is_selectable_not_outlined_paths():
    # svg.fonttype must stay "none" so labels/sequences are real <text> elements
    # (selectable, copyable) rather than vector-outlined glyphs.
    assert matplotlib.rcParams["svg.fonttype"] == "none"


def test_hit_stats_text_includes_containment_jaccard_enrichment_pvalue():
    # Both scopes are shown: either one can be why the hit was reported at all.
    hit = vh._build_hit("tgt", [_row(containment=0.5, jaccard=0.123, region_enrichment=4.5,
                                     region_poisson_pvalue=1e-06, query_poisson_pvalue=0.99)])
    stats = vh.GenePlot._hit_stats_text(hit)
    assert "containment=0.50" in stats
    assert "jaccard=0.123" in stats
    assert "region enrich=4.50" in stats
    assert "region p=1e-06" in stats
    assert "query p=0.99" in stats


def test_plot_gene_svg_contains_stats_line(tmp_path):
    hits = vh.merge_regions_by_target(
        [_row(target_name="TGT_A", region_start=10, region_end=40, region_length=30,
              containment=0.5, jaccard=0.123, region_enrichment=4.5,
              region_poisson_pvalue=1e-06)],
        gap_merge=10,
    )
    svg_path = tmp_path / "gene.hits.svg"
    vh.plot_gene(CED9_NAME, CED9_LEN, hits, [str(svg_path)])
    svg_text = svg_path.read_text()
    assert "jaccard=0.123" in svg_text
    assert "region enrich=4.50" in svg_text
    assert "region p=1e-06" in svg_text


# --- main() end-to-end CLI ----------------------------------------------------

def _write_csv(path, rows):
    # Derived from _row() rather than listed again, so there is exactly one place
    # in this file that has to track the CSV schema.
    fieldnames = list(_row())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            full_row = {k: row.get(k, "") for k in fieldnames}
            writer.writerow(full_row)


def test_main_end_to_end_writes_one_png_per_query(tmp_path, capsys):
    csv_path = tmp_path / "results.csv"
    _write_csv(csv_path, [
        _row(query_name=SEPT4_NAME, target_name="tgt1", region_start=0, region_end=10, region_length=10),
        _row(query_name=B2L11_NAME, target_name="tgt2", region_start=0, region_end=10, region_length=10),
    ])

    out_dir = tmp_path / "out"
    sys.argv = ["visualize_hits.py", "--csv", str(csv_path), "--query-fasta", MULTI_FASTA,
                "--output-dir", str(out_dir)]
    vh.main()

    assert (out_dir / "SEPT4_HUMAN.hits.png").exists()
    assert (out_dir / "SEPT4_HUMAN.hits.svg").exists()
    assert (out_dir / "B2L11_HUMAN.hits.png").exists()
    assert (out_dir / "B2L11_HUMAN.hits.svg").exists()
    out = capsys.readouterr().out
    assert "Wrote" in out and "SEPT4_HUMAN.hits.png" in out
    assert "1 targets, 1 hits, 1 regions" in out


def test_main_skips_query_missing_from_fasta(tmp_path, capsys):
    # BCL2_NAME is a real header, just not one that's in CED9_FASTA.
    csv_path = tmp_path / "results.csv"
    _write_csv(csv_path, [_row(query_name=BCL2_NAME, region_start=0, region_end=10)])

    out_dir = tmp_path / "out"
    sys.argv = ["visualize_hits.py", "--csv", str(csv_path), "--query-fasta", CED9_FASTA,
                "--output-dir", str(out_dir)]
    vh.main()

    assert list(out_dir.glob("*.png")) == []
    out = capsys.readouterr().out
    assert f"Skipping '{BCL2_NAME}': not found" in out


def test_main_min_containment_filters_low_confidence_rows(tmp_path):
    csv_path = tmp_path / "results.csv"
    _write_csv(csv_path, [_row(query_name=CED9_NAME, containment=0.1, region_start=0, region_end=10)])

    out_dir = tmp_path / "out"
    sys.argv = ["visualize_hits.py", "--csv", str(csv_path), "--query-fasta", CED9_FASTA,
                "--output-dir", str(out_dir), "--min-containment", "0.5"]
    vh.main()

    assert list(out_dir.glob("*.png")) == []


def test_main_max_hits_caps_by_distinct_target_not_by_fragment_count(tmp_path, capsys):
    # "frag" has 3 widely-spaced (unmergeable) hits at high containment; "other" has
    # a single hit at lower containment. --max-hits 1 must keep only the best distinct
    # target ("frag", all 3 of its fragments) rather than being consumed by fragment count.
    csv_path = tmp_path / "results.csv"
    _write_csv(csv_path, [
        _row(query_name=CED9_NAME, target_name="frag", containment=0.9, region_start=0, region_end=10),
        _row(query_name=CED9_NAME, target_name="frag", containment=0.9, region_start=100, region_end=110),
        _row(query_name=CED9_NAME, target_name="frag", containment=0.9, region_start=200, region_end=210),
        _row(query_name=CED9_NAME, target_name="other", containment=0.5, region_start=50, region_end=60),
    ])

    out_dir = tmp_path / "out"
    sys.argv = ["visualize_hits.py", "--csv", str(csv_path), "--query-fasta", CED9_FASTA,
                "--output-dir", str(out_dir), "--max-hits", "1"]
    vh.main()

    out = capsys.readouterr().out
    assert "1 targets, 3 hits, 4 regions" in out


def test_main_corrects_pvalue_across_all_targets_not_just_displayed_ones(tmp_path):
    # "hidden" is a real target for this query but gets filtered out by
    # --min-containment; its p-value must still count toward the correction
    # denominator for the targets that *are* displayed, or the q-value would
    # understate how many comparisons were actually made.
    csv_path = tmp_path / "results.csv"
    _write_csv(csv_path, [
        _row(query_name=CED9_NAME, target_name="shown", containment=0.9,
             region_poisson_pvalue=0.01, region_start=0, region_end=10),
        _row(query_name=CED9_NAME, target_name="hidden", containment=0.01,
             region_poisson_pvalue=0.02, region_start=50, region_end=60),
    ])

    out_dir = tmp_path / "out"
    sys.argv = ["visualize_hits.py", "--csv", str(csv_path), "--query-fasta", CED9_FASTA,
                "--output-dir", str(out_dir), "--min-containment", "0.5"]
    vh.main()

    svg_text = (out_dir / f"{vh.safe_filename(CED9_NAME)}.hits.svg").read_text()
    expected_q = vh.benjamini_hochberg({"shown": 0.01, "hidden": 0.02})["shown"]
    assert f"region q={expected_q:.2g}" in svg_text


def test_main_query_name_filters_to_single_gene(tmp_path):
    csv_path = tmp_path / "results.csv"
    _write_csv(csv_path, [
        _row(query_name=SEPT4_NAME, region_start=0, region_end=10),
        _row(query_name=B2L11_NAME, region_start=0, region_end=10),
    ])

    out_dir = tmp_path / "out"
    sys.argv = ["visualize_hits.py", "--csv", str(csv_path), "--query-fasta", MULTI_FASTA,
                "--output-dir", str(out_dir), "--query-name", SEPT4_NAME]
    vh.main()

    assert (out_dir / "SEPT4_HUMAN.hits.png").exists()
    assert not (out_dir / "B2L11_HUMAN.hits.png").exists()


def test_main_query_name_accepts_short_gene_symbol(tmp_path):
    csv_path = tmp_path / "results.csv"
    _write_csv(csv_path, [
        _row(query_name=SEPT4_NAME, region_start=0, region_end=10),
        _row(query_name=B2L11_NAME, region_start=0, region_end=10),
    ])

    out_dir = tmp_path / "out"
    sys.argv = ["visualize_hits.py", "--csv", str(csv_path), "--query-fasta", MULTI_FASTA,
                "--output-dir", str(out_dir), "--query-name", "SEPT4"]
    vh.main()

    assert (out_dir / "SEPT4_HUMAN.hits.png").exists()
    assert not (out_dir / "B2L11_HUMAN.hits.png").exists()


def test_main_query_name_no_match_prints_available_genes(tmp_path, capsys):
    csv_path = tmp_path / "results.csv"
    _write_csv(csv_path, [_row(query_name=SEPT4_NAME, region_start=0, region_end=10)])

    out_dir = tmp_path / "out"
    sys.argv = ["visualize_hits.py", "--csv", str(csv_path), "--query-fasta", MULTI_FASTA,
                "--output-dir", str(out_dir), "--query-name", "NOT_A_GENE"]
    vh.main()

    assert list(out_dir.glob("*.png")) == []
    out = capsys.readouterr().out
    assert "No query matching 'NOT_A_GENE' found" in out
    assert "SEPT4_HUMAN" in out
