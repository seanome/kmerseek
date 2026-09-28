#!/usr/bin/env python3
"""Render one interactive HTML report per query from `kmerseek search` output.

The query is the shared axis. It is drawn once at the top as a line with its domains
(from --domains), under a histogram of how many database entries have a run over each
residue: light for any run, dark for a run with --solid-identical or more identical
residues. That histogram is the noise map: a low-complexity stretch such as the BCL-2
loop is covered by a third of unrelated proteins, so a run there is discounted at a
glance and a run in BH1 is not.

Below it, one row per protein with the numbers in the row (length, runs, longest run,
identical residues in it, shared k-mers, the ranking statistic, and with --structures
the TM-score) and every run drawn as a bar at its query coordinates, shaded by its
share of identical residues; overlapping runs stack in lanes. Database entries of one
gene (UniProt GN= and OS=) fold into one row, so a family search is not a list of
TrEMBL copies of the query. Clicking a row opens the pair view underneath it: the dot
plot with protein tracks and one alignment block per run.

Rows are ordered by `region_evalue` when the CSV has it (kmerseek >= 0.5) and otherwise
by the Benjamini-Hochberg corrected region tail probability. Sorting by identical
residues or run length instead puts composition-driven hits (p53, POU4F1) among family
members, which is why the ranking statistic is the default. The page filters rows by
name, identical residues, the statistic, Swiss-Prot status, fragments and which query
domain a run falls in, and downloads the rows and runs as TSV, the runs as FASTA, the
overview as SVG or PNG, and its data as JSON.

The page itself is `kmerseek_hits_template.html` next to this script (see hits_page.py);
the template's own script draws everything from the data.

The pair view needs every shared k-mer, which the CSV does not carry, so this script runs
`kmerseek pair` once per row on sequences taken from the two FASTA files.

Usage:
    kmerseek search -q queries.fasta -t targets.rocksdb -o results.csv --alphabet hp --ksize 12
    python visualize_search.py --csv results.csv --query-fasta queries.fasta \\
        --target-fasta targets.fasta.gz --output-dir report/ --domains pfam_domains.tsv
"""

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pair_model import build_model, count_agreement, domains_for, load_domains
from structure_alignment import align_pair, find_aligner
from visualize_hits import (
    _target_best_rows,
    _target_tail_probabilities,
    benjamini_hochberg,
    iter_query_rows,
    load_query_names,
    load_rows_for_queries,
    resolve_query_names,
    safe_filename,
    scan_csv,
    short_label,
)
from hits_page import render_page

# -- sequences --------------------------------------------------------------------------


def read_fasta(path):
    """{header after '>': sequence}, gzip or plain."""
    opener = gzip.open if path.endswith(".gz") else open
    records, name, parts = {}, None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts)
                name, parts = line[1:], []
            else:
                parts.append(line.strip())
    if name is not None:
        records[name] = "".join(parts)
    return records


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n{seq}\n")


def has_pair(binary):
    """Whether this binary knows `kmerseek pair` (older releases do not)."""
    try:
        return subprocess.run([binary, "pair", "--help"], capture_output=True).returncode == 0
    except OSError:
        return False


def find_kmerseek(explicit):
    """The kmerseek binary: --kmerseek, then PATH, then this checkout's release or debug
    build, skipping any that predates `kmerseek pair`."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    candidates = [explicit, shutil.which("kmerseek"), os.path.join(root, "target", "release", "kmerseek"), os.path.join(root, "target", "debug", "kmerseek")]
    for c in candidates:
        if c and os.path.exists(c) and has_pair(c):
            return c
    raise SystemExit("no kmerseek binary with the `pair` subcommand found; pass --kmerseek or build with `cargo build --release`")


def run_pair(kmerseek, query_fasta, target_fasta, target_name, ksize, moltype):
    """`kmerseek pair` for one hit, as the parsed JSON."""
    cmd = [kmerseek, "pair", "--query", query_fasta, "--target", target_fasta, "--target-name", target_name, "--ksize", str(ksize), "--alphabet", moltype]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        raise SystemExit(f"kmerseek pair failed for {target_name}:\n{out.stderr.strip()}")
    return json.loads(out.stdout)


# -- headers and ranking --------------------------------------------------------------------


def header_field(header, key):
    """`GN=BCL2` style field of a UniProt header, or None."""
    m = re.search(rf"(?:^|\s){key}=(.+?)(?=\s\w\w?=|$)", header)
    return m.group(1) if m else None


def description(header):
    """The FASTA header after its first token, with UniProt's OS=... tail dropped."""
    rest = header.split(" ", 1)[1] if " " in header else ""
    return rest.split(" OS=")[0]


def protein_key(header):
    """What folds database entries into one protein row: gene and organism when the
    header has them (gene lower-cased with punctuation dropped, since TrEMBL writes
    bcl-2 for BCL2), else the header itself."""
    gene, organism = header_field(header, "GN"), header_field(header, "OS")
    if not gene:
        return (header, None)
    return (re.sub(r"[^a-z0-9]", "", gene.lower()), organism)


def ranking(rows):
    """{target_name: (statistic, is_evalue)} for ordering rows: region_evalue when the CSV
    carries it, else the BH q-value of the best region's tail probability."""
    best = _target_best_rows(rows)
    if "region_evalue" in rows[0]:
        return {t: float(r["region_evalue"]) for t, r in best.items()}, True
    return benjamini_hochberg(_target_tail_probabilities(rows)), False


def fold_entries(rows):
    """{protein key: [target names]} for every target in the CSV."""
    groups = defaultdict(list)
    for name in _target_best_rows(rows):
        groups[protein_key(name)].append(name)
    return groups


def rank_proteins(rows, max_rows):
    """[(representative target name, other entry names, statistic)] one per protein,
    best statistic first, capped at max_rows."""
    stat, _ = ranking(rows)
    best = _target_best_rows(rows)
    # Best statistic first, then the higher region score; on a full tie the reviewed (sp|)
    # entry represents the protein.
    order = lambda n: (stat[n], -float(best[n]["region_poisson_score"]), not n.startswith("sp|"))
    proteins = []
    for names in fold_entries(rows).values():
        names.sort(key=order)
        proteins.append((names[0], names[1:], stat[names[0]]))
    proteins.sort(key=lambda p: order(p[0]))
    return proteins[:max_rows]


# -- histogram over the query ----------------------------------------------------------------


def region_identical(row):
    return count_agreement(row["region_subseq"], row["target_subseq"])


def coverage(rows, query_length, ksize, solid_identical):
    """Per query residue (0-based), how many database entries have a run over it, and how
    many have one with at least solid_identical identical residues. From the CSV alone, so
    it counts every entry the search reported, not only the rows shown."""
    any_run, solid = [0] * query_length, [0] * query_length
    per_target = defaultdict(lambda: (set(), set()))
    for row in rows:
        if int(row["region_length"]) <= ksize:
            continue
        span = range(int(row["region_start"]), int(row["region_end"]))
        covered, covered_solid = per_target[row["target_name"]]
        covered.update(span)
        if region_identical(row) >= solid_identical:
            covered_solid.update(span)
    for covered, covered_solid in per_target.values():
        for i in covered:
            any_run[i] += 1
        for i in covered_solid:
            solid[i] += 1
    return {"any": any_run, "solid": solid}


# -- rows --------------------------------------------------------------------------------------


def run_identical(block):
    return block["identical"]


def run_bar(block):
    """What a row's bar and its tooltip need, 0-based half-open as in the model."""
    bar = {k: block[k] for k in ("number", "query_start", "query_end", "target_start", "target_end", "length", "polar")}
    bar["identical"] = run_identical(block)
    if "structure_offset" in block:
        bar["structure_offset"] = block["structure_offset"]
    return bar


def protein_row(rank, target_name, others, stat, row, model):
    best = model["runs"][0] if model["runs"] else None
    return {
        "rank": rank,
        "label": model["target"]["label"],
        "target_name": target_name,
        "entry": target_name.split()[0],
        "gene": header_field(target_name, "GN"),
        "organism": header_field(target_name, "OS"),
        "description": description(target_name),
        "n_entries": 1 + len(others),
        "other_entries": [o.split()[0] for o in others],
        "length": model["target"]["length"],
        "n_runs": len(model["runs"]),
        "best_length": best["length"] if best else 0,
        "best_identical": run_identical(best) if best else 0,
        "tm_score": model["structure"]["tm_score_query"] if model.get("structure") else None,
        "n_shared": int(row["n_intersecting_hashes"]),
        "stat": stat,
        "runs": [run_bar(b) for b in model["runs"]],
        "model": model,
    }


class SearchReport:
    """Builds the per-query report, running `kmerseek pair` for each protein row."""

    def __init__(self, args, kmerseek, domain_rows):
        self.args = args
        self.kmerseek = kmerseek
        self.domain_rows = domain_rows
        self.queries = read_fasta(args.query_fasta)
        self.targets = read_fasta(args.target_fasta)
        self.aligner = find_aligner(args.aligner) if args.structures else None
        if args.structures and self.aligner is None:
            print("no USalign or TMalign found; skipping the superpositions", file=sys.stderr)

    def structure(self, query_name, target_name):
        if not (self.args.structures and self.aligner):
            return None
        return align_pair(self.aligner, self.args.structures, query_name, target_name)

    def rows_for(self, query_name, rows, workdir):
        ranked = rank_proteins(rows, self.args.max_rows)
        best = _target_best_rows(rows)
        ksize, moltype = int(rows[0]["ksize"]), rows[0]["moltype"]
        query_fasta = os.path.join(workdir, "query.fasta")
        write_fasta(query_fasta, {query_name: self.queries[query_name]})
        target_fasta = os.path.join(workdir, "targets.fasta")
        write_fasta(target_fasta, {t: self.targets[t] for t, _, _ in ranked})
        out = []
        for rank, (target_name, others, stat) in enumerate(ranked, start=1):
            pair = run_pair(self.kmerseek, query_fasta, target_fasta, target_name.split()[0], ksize, moltype)
            model = build_model(
                pair,
                self.domain_rows,
                flank=self.args.flank,
                structure=self.structure(query_name, target_name),
            )
            out.append(protein_row(rank, target_name, others, stat, best[target_name], model))
        return out

    def query_track(self, query_name):
        return {
            "label": short_label(query_name),
            "name": query_name,
            "length": len(self.queries[query_name]),
            "domains": domains_for(self.domain_rows, query_name),
        }

    def report(self, query_name, rows):
        with tempfile.TemporaryDirectory() as workdir:
            protein_rows = self.rows_for(query_name, rows, workdir)
        ksize = int(rows[0]["ksize"])
        _, is_evalue = ranking(rows)
        query = self.query_track(query_name)
        return {
            "query": query,
            "ksize": ksize,
            "moltype": rows[0]["moltype"],
            "alphabet": protein_rows[0]["model"]["alphabet"] if protein_rows else rows[0]["moltype"],
            "classes": protein_rows[0]["model"]["classes"] if protein_rows else [],
            "stat_name": "E-value" if is_evalue else "q-value",
            "stat_note": "Karlin-Altschul E-value of the best region" if is_evalue else "Benjamini-Hochberg corrected tail probability of the best region",
            "n_entries": len(_target_best_rows(rows)),
            "n_proteins": len(fold_entries(rows)),
            "solid_identical": self.args.solid_identical,
            "structures": bool(self.args.structures and self.aligner),
            "max_runs_shown": self.args.max_runs_shown,
            "coverage": coverage(rows, query["length"], ksize, self.args.solid_identical),
            "rows": protein_rows,
        }

    def render(self, query_name, rows):
        return render_report(self.report(query_name, rows))


# -- HTML ----------------------------------------------------------------------------------


def render_report(report):
    return render_page(f"{report['query']['label']} kmerseek hits", report)


# -- CLI -----------------------------------------------------------------------------------


def _build_arg_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv", required=True, help="CSV written by `kmerseek search -o`")
    p.add_argument("--query-fasta", required=True, help="the FASTA the search was run with")
    p.add_argument("--target-fasta", required=True, help="the FASTA the target index was built from (gzip or plain)")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--query-name", help="one query's header, or its first token; default every query in the CSV")
    p.add_argument("--domains", nargs="*", default=[], metavar="TABLE", help="Pfam-style domain tables for queries and targets")
    p.add_argument("--max-rows", type=int, default=100, help="protein rows per query, best ranking statistic first (default 100)")
    p.add_argument("--max-runs-shown", type=int, default=10, help="alignments per opened row, longest first (default 10)")
    p.add_argument("--solid-identical", type=int, default=5, help="identical residues a run needs to count in the histogram's dark area (default 5)")
    p.add_argument("--flank", type=int, default=0, help="residues shown either side of each run in the alignments")
    p.add_argument("--structures", metavar="DIR", help="directory of AlphaFold or PDB files; with USalign or TM-align, each row gets a TM-score and its dot plot their residue pairs")
    p.add_argument("--aligner", help="USalign or TMalign binary (default: found on PATH)")
    p.add_argument("--kmerseek", help="path to the kmerseek binary (default: PATH, then target/release, target/debug)")
    return p


def main():
    args = _build_arg_parser().parse_args()
    kmerseek = find_kmerseek(args.kmerseek)
    lazy = scan_csv(args.csv)
    names = resolve_query_names(args.query_name, load_query_names(lazy))
    report = SearchReport(args, kmerseek, load_domains(args.domains))
    os.makedirs(args.output_dir, exist_ok=True)
    for query_name, rows in iter_query_rows(load_rows_for_queries(lazy, names)):
        path = os.path.join(args.output_dir, f"{safe_filename(query_name)}.kmerseek_hits.html")
        with open(path, "w") as fh:
            fh.write(report.render(query_name, rows))
        print(path)


if __name__ == "__main__":
    main()
