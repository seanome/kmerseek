"""Tests for visualize_hits.py.

Run with: /Users/olga/anaconda3/envs/kmerseek-dev/bin/python3 -m pytest scripts/test_visualize_hits.py -v
"""

import csv
import os
import sys

import matplotlib

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

def _row(target_name="tgt", query_start=0, query_end=10, region_length=10,
         containment=0.5, moltype="hp", query_subseq="MKVLLLKKKK",
         moltype_seq="hpphhhpppp", target_subseq="TTTTTTTTTT", query_name=CED9_NAME):
    return {
        "query_name": query_name,
        "target_name": target_name,
        "query_start": str(query_start),
        "query_end": str(query_end),
        "region_length": str(region_length),
        "containment": str(containment),
        "moltype": moltype,
        "query_subseq": query_subseq,
        "moltype_seq": moltype_seq,
        "target_subseq": target_subseq,
    }


def test_merge_regions_within_gap_becomes_one_hit():
    rows = [
        _row(query_start=0, query_end=10, region_length=10),
        _row(query_start=15, query_end=25, region_length=10),  # gap of 5
    ]
    hits = vh.merge_regions_by_target(rows, gap_merge=10)
    assert len(hits) == 1
    assert hits[0]["start"] == 0
    assert hits[0]["end"] == 25
    assert hits[0]["n_regions"] == 2


def test_merge_regions_beyond_gap_stays_separate():
    rows = [
        _row(query_start=0, query_end=10, region_length=10),
        _row(query_start=25, query_end=35, region_length=10),  # gap of 15
    ]
    hits = vh.merge_regions_by_target(rows, gap_merge=10)
    assert len(hits) == 2
    assert (hits[0]["start"], hits[0]["end"]) == (0, 10)
    assert (hits[1]["start"], hits[1]["end"]) == (25, 35)


def test_merge_regions_groups_by_target_independently():
    rows = [
        _row(target_name="A", query_start=0, query_end=10, region_length=10),
        _row(target_name="B", query_start=5, query_end=15, region_length=10),
    ]
    hits = vh.merge_regions_by_target(rows, gap_merge=10)
    assert {h["target_name"] for h in hits} == {"A", "B"}
    assert len(hits) == 2


def test_build_hit_region_rows_keeps_every_region_sorted_by_start():
    # Every region must stay visible (not just one "representative" pick) so a
    # multi-region hit like BAK_HUMAN's two separate matches both get shown.
    second_region = _row(query_start=20, query_end=25, region_length=5,
                          query_subseq="MKVLL", moltype_seq="hpphh",
                          target_subseq="AAAAA", containment=0.9)
    first_region = _row(query_start=0, query_end=20, region_length=20,
                         query_subseq="MKVLLLKKKKMKVLLLKKKK",
                         moltype_seq="hpphhhppppphhhhhppppp"[:20],
                         target_subseq="TTTTTTTTTTTTTTTTTTTT", containment=0.3)
    hit = vh._build_hit("tgt", [second_region, first_region])
    assert [row["query_subseq"] for row in hit["region_rows"]] == [
        "MKVLLLKKKKMKVLLLKKKK", "MKVLL",
    ]


def test_build_hit_containment_is_max_across_cluster():
    rows = [_row(containment=0.2), _row(containment=0.7), _row(containment=0.4)]
    hit = vh._build_hit("tgt", rows)
    assert hit["containment"] == 0.7


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
        _row(query_start=0, query_end=10, region_length=10),
        _row(query_start=5, query_end=15, region_length=10),
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
        [_row(target_name="TGT_A", query_start=10, query_end=40, region_length=30)],
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


# --- main() end-to-end CLI ----------------------------------------------------

def _write_csv(path, rows):
    fieldnames = ["query_name", "target_name", "containment", "query_start", "query_end",
                  "query_subseq", "target_start", "target_end", "target_subseq",
                  "moltype_seq", "moltype", "region_length"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            full_row = {k: row.get(k, "") for k in fieldnames}
            writer.writerow(full_row)


def test_main_end_to_end_writes_one_png_per_query(tmp_path, capsys):
    csv_path = tmp_path / "results.csv"
    _write_csv(csv_path, [
        _row(query_name=SEPT4_NAME, target_name="tgt1", query_start=0, query_end=10, region_length=10),
        _row(query_name=B2L11_NAME, target_name="tgt2", query_start=0, query_end=10, region_length=10),
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
    _write_csv(csv_path, [_row(query_name=BCL2_NAME, query_start=0, query_end=10)])

    out_dir = tmp_path / "out"
    sys.argv = ["visualize_hits.py", "--csv", str(csv_path), "--query-fasta", CED9_FASTA,
                "--output-dir", str(out_dir)]
    vh.main()

    assert list(out_dir.glob("*.png")) == []
    out = capsys.readouterr().out
    assert f"Skipping '{BCL2_NAME}': not found" in out


def test_main_min_containment_filters_low_confidence_rows(tmp_path):
    csv_path = tmp_path / "results.csv"
    _write_csv(csv_path, [_row(query_name=CED9_NAME, containment=0.1, query_start=0, query_end=10)])

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
        _row(query_name=CED9_NAME, target_name="frag", containment=0.9, query_start=0, query_end=10),
        _row(query_name=CED9_NAME, target_name="frag", containment=0.9, query_start=100, query_end=110),
        _row(query_name=CED9_NAME, target_name="frag", containment=0.9, query_start=200, query_end=210),
        _row(query_name=CED9_NAME, target_name="other", containment=0.5, query_start=50, query_end=60),
    ])

    out_dir = tmp_path / "out"
    sys.argv = ["visualize_hits.py", "--csv", str(csv_path), "--query-fasta", CED9_FASTA,
                "--output-dir", str(out_dir), "--max-hits", "1"]
    vh.main()

    out = capsys.readouterr().out
    assert "1 targets, 3 hits, 4 regions" in out


def test_main_query_name_filters_to_single_gene(tmp_path):
    csv_path = tmp_path / "results.csv"
    _write_csv(csv_path, [
        _row(query_name=SEPT4_NAME, query_start=0, query_end=10),
        _row(query_name=B2L11_NAME, query_start=0, query_end=10),
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
        _row(query_name=SEPT4_NAME, query_start=0, query_end=10),
        _row(query_name=B2L11_NAME, query_start=0, query_end=10),
    ])

    out_dir = tmp_path / "out"
    sys.argv = ["visualize_hits.py", "--csv", str(csv_path), "--query-fasta", MULTI_FASTA,
                "--output-dir", str(out_dir), "--query-name", "SEPT4"]
    vh.main()

    assert (out_dir / "SEPT4_HUMAN.hits.png").exists()
    assert not (out_dir / "B2L11_HUMAN.hits.png").exists()


def test_main_query_name_no_match_prints_available_genes(tmp_path, capsys):
    csv_path = tmp_path / "results.csv"
    _write_csv(csv_path, [_row(query_name=SEPT4_NAME, query_start=0, query_end=10)])

    out_dir = tmp_path / "out"
    sys.argv = ["visualize_hits.py", "--csv", str(csv_path), "--query-fasta", MULTI_FASTA,
                "--output-dir", str(out_dir), "--query-name", "NOT_A_GENE"]
    vh.main()

    assert list(out_dir.glob("*.png")) == []
    out = capsys.readouterr().out
    assert "No query matching 'NOT_A_GENE' found" in out
    assert "SEPT4_HUMAN" in out
