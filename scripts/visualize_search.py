#!/usr/bin/env python3
"""Render one interactive HTML report per query from `kmerseek search` output, laid out
like a Foldseek results page.

One row per target hit, ordered by Benjamini-Hochberg corrected region tail probability
(q-value), with the numbers in the row: runs of consecutive shared k-mers, shared k-mers,
containment, whole-query p-value, best region score, q-value, and a bar showing where
the hit's runs land on the query. Clicking a row opens the pair view underneath it: the
dot plot with protein tracks and every run's alignment, the same panel
visualize_pair.py draws.

The pair view needs every shared k-mer, which the search CSV does not carry, so this
script runs `kmerseek pair` once per hit on the query and target sequences taken from
the two FASTA files. --domains takes the same Pfam-style tables as visualize_pair.py and
labels both proteins.

Usage:
    kmerseek search -q queries.fasta -t targets.rocksdb -o results.csv --alphabet hp --ksize 12
    python visualize_search.py --csv results.csv --query-fasta queries.fasta \\
        --target-fasta targets.fasta.gz --output-dir report/ --domains pfam_domains.tsv
"""

import argparse
import gzip
import html
import json
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pair_model import build_model, domains_for, load_domains
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
from visualize_pair_html import PAIR_CSS, PAIR_JS, embed_json

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


# -- hits ----------------------------------------------------------------------------------


def rank_targets(rows, max_hits):
    """[(target_name, best row, q-value)] ordered by q-value then best region score,
    capped at max_hits distinct targets."""
    best = _target_best_rows(rows)
    q = benjamini_hochberg(_target_tail_probabilities(rows))
    ranked = sorted(best, key=lambda t: (q[t], -float(best[t]["region_poisson_score"])))
    return [(t, best[t], q[t]) for t in ranked[:max_hits]]


def description(header):
    """The FASTA header after its first token, with UniProt's OS=... tail dropped."""
    rest = header.split(" ", 1)[1] if " " in header else ""
    return rest.split(" OS=")[0]


def hit_entry(rank, target_name, row, q_value, model):
    return {
        "rank": rank,
        "target": model["target"]["label"],
        "target_name": target_name,
        "description": description(target_name),
        "n_runs": len(model["runs"]),
        "n_shared": int(row["n_intersecting_hashes"]),
        "containment": float(row["containment"]),
        "query_pvalue": float(row["query_poisson_pvalue"]),
        "region_score": float(row["region_poisson_score"]),
        "q_value": q_value,
        "model": model,
    }


class SearchReport:
    """Builds the per-query hit list, running `kmerseek pair` for each hit."""

    def __init__(self, args, kmerseek, domain_rows):
        self.args = args
        self.kmerseek = kmerseek
        self.domain_rows = domain_rows
        self.queries = read_fasta(args.query_fasta)
        self.targets = read_fasta(args.target_fasta)

    def hits_for(self, query_name, rows, workdir):
        ranked = rank_targets(rows, self.args.max_hits)
        ksize, moltype = int(rows[0]["ksize"]), rows[0]["moltype"]
        query_fasta = os.path.join(workdir, "query.fasta")
        write_fasta(query_fasta, {query_name: self.queries[query_name]})
        target_fasta = os.path.join(workdir, "targets.fasta")
        write_fasta(target_fasta, {t: self.targets[t] for t, _, _ in ranked})
        hits = []
        for rank, (target_name, row, q_value) in enumerate(ranked, start=1):
            pair = run_pair(self.kmerseek, query_fasta, target_fasta, target_name.split()[0], ksize, moltype)
            model = build_model(pair, self.domain_rows, flank=self.args.flank)
            hits.append(hit_entry(rank, target_name, row, q_value, model))
        return hits

    def query_track(self, query_name):
        return {
            "label": short_label(query_name),
            "name": query_name,
            "length": len(self.queries[query_name]),
            "domains": domains_for(self.domain_rows, query_name),
        }

    def render(self, query_name, rows):
        with tempfile.TemporaryDirectory() as workdir:
            hits = self.hits_for(query_name, rows, workdir)
        report = {
            "query": self.query_track(query_name),
            "ksize": int(rows[0]["ksize"]),
            "moltype": rows[0]["moltype"],
            "alphabet": hits[0]["model"]["alphabet"] if hits else rows[0]["moltype"],
            "n_targets_tested": len(_target_best_rows(rows)),
            "hits": hits,
        }
        return render_report(report)


# -- HTML ----------------------------------------------------------------------------------

REPORT_CSS = r"""
  h1 { font-size: 16px; font-weight: 600; margin: 0 0 2px; }
  .meta { color: var(--secondary); margin: 0 0 12px; }
  .qtrack { margin: 4px 0 14px; }
  table { border-collapse: collapse; width: 100%; font-size: 13px; }
  th, td { text-align: left; padding: 6px 10px; border-bottom: 1px solid #e3e2dc; vertical-align: middle; white-space: nowrap; }
  th { color: var(--secondary); font-weight: 600; cursor: pointer; user-select: none; }
  th.num, td.num { text-align: right; font-variant-numeric: tabular-nums; }
  th.sorted:after { content: " \25BE"; }
  th.sorted.asc:after { content: " \25B4"; }
  td.desc { max-width: 260px; overflow: hidden; text-overflow: ellipsis; color: var(--secondary); }
  tr.hit { cursor: pointer; }
  tr.hit:hover { background: #f1f0eb; }
  tr.hit.open { background: #ecebe4; }
  tr.detail td { padding: 12px 10px 18px; white-space: normal; background: #fafaf7; }
  .close { float: right; font: inherit; border: 1px solid var(--edge); background: #fff; border-radius: 4px; padding: 2px 8px; cursor: pointer; }
  .bar svg { display: block; }
  .scroll { overflow-x: auto; }
"""

REPORT_JS = r"""
const R = __REPORT__;
const BAR_W = 220, BAR_H = 22;

function positionBar(track, model) {
  // The query as a line with its domain boxes; each run of the hit as a dark segment on it.
  const scale = BAR_W / track.length;
  const svg = svgEl("svg", { width: BAR_W + 44, height: BAR_H });
  const x = p => 22 + (p - 1) * scale;
  svg.appendChild(svgEl("line", { x1: x(1), x2: x(track.length), y1: 11, y2: 11, stroke: "var(--domain-edge)" }));
  for (const d of track.domains)
    svg.appendChild(svgEl("rect", { x: x(d.start), y: 6, width: (d.end - d.start + 1) * scale, height: 10, fill: "var(--domain)", stroke: "var(--domain-edge)" }));
  for (const s of model.singles)
    svg.appendChild(svgEl("circle", { cx: x(s.query_pos + (model.ksize + 1) / 2), cy: 11, r: 1.6, fill: "var(--single)" }));
  for (const b of model.runs)
    svg.appendChild(svgEl("line", { x1: x(b.query_start + 1), x2: x(b.query_end), y1: 11, y2: 11, stroke: "var(--run)", "stroke-width": 4, "stroke-linecap": "round" }));
  svg.appendChild(svgEl("text", { x: 18, y: 15, "text-anchor": "end", "font-size": 9 }, 1));
  svg.appendChild(svgEl("text", { x: x(track.length) + 4, y: 15, "font-size": 9 }, track.length));
  return svg;
}

function queryTrack(track) {
  // The query alone, wider, with domain names, above the table.
  const W = 520, scale = W / track.length;
  const svg = svgEl("svg", { width: W + 60, height: 40 });
  const x = p => 30 + (p - 1) * scale;
  svg.appendChild(svgEl("line", { x1: x(1), x2: x(track.length), y1: 26, y2: 26, stroke: "var(--domain-edge)" }));
  track.domains.forEach((d, i) => {
    svg.appendChild(svgEl("rect", { x: x(d.start), y: 20, width: (d.end - d.start + 1) * scale, height: 12, fill: "var(--domain)", stroke: "var(--domain-edge)" }));
    svg.appendChild(svgEl("text", { x: x((d.start + d.end) / 2), y: 14 - (labelsCollide(track, scale) && i % 2 ? 11 : 0), "text-anchor": "middle" }, d.name));
  });
  svg.appendChild(svgEl("text", { x: 26, y: 30, "text-anchor": "end", "font-size": 9 }, 1));
  svg.appendChild(svgEl("text", { x: x(track.length) + 4, y: 30, "font-size": 9 }, track.length));
  return svg;
}

const fmtP = p => p === 0 ? "0" : p < 1e-3 ? p.toExponential(2) : p.toFixed(3);
const COLUMNS = [
  ["rank", "#", "num"], ["target", "Target", ""], ["description", "Description", "desc"],
  ["n_runs", "Runs", "num"], ["n_shared", `Shared ${R.ksize}-mers`, "num"], ["containment", "Containment", "num"],
  ["query_pvalue", "Query p-value", "num"], ["region_score", "Best region score", "num"], ["q_value", "q-value", "num"],
  ["bar", "Position in query", ""],
];
let sortKey = "rank", sortAsc = true, openRank = null;

function cellText(h, key) {
  if (key === "containment") return h.containment.toFixed(3);
  if (key === "query_pvalue" || key === "q_value") return fmtP(h[key]);
  if (key === "region_score") return h.region_score.toFixed(2);
  return String(h[key]);
}

function renderTable() {
  const table = document.getElementById("hits");
  table.innerHTML = "";
  const head = el("tr");
  for (const [key, label, cls] of COLUMNS) {
    const th = el("th", cls + (key === sortKey ? " sorted" + (sortAsc ? " asc" : "") : ""), label);
    if (key !== "bar") th.onclick = () => { sortAsc = key === sortKey ? !sortAsc : key === "target" || key === "description"; sortKey = key; renderTable(); };
    head.appendChild(th);
  }
  table.appendChild(head);
  const hits = [...R.hits].sort((a, b) => (a[sortKey] < b[sortKey] ? -1 : a[sortKey] > b[sortKey] ? 1 : 0) * (sortAsc ? 1 : -1));
  for (const h of hits) {
    const tr = el("tr", "hit" + (h.rank === openRank ? " open" : ""));
    for (const [key, , cls] of COLUMNS) {
      const td = el("td", cls);
      if (key === "bar") { td.className = "bar"; td.appendChild(positionBar(R.query, h.model)); }
      else td.textContent = cellText(h, key);
      if (key === "description") td.title = h.target_name;
      tr.appendChild(td);
    }
    tr.onclick = () => { openRank = openRank === h.rank ? null : h.rank; renderTable(); };
    table.appendChild(tr);
    if (h.rank === openRank) table.appendChild(detailRow(h));
  }
}

function detailRow(h) {
  const tr = el("tr", "detail"), td = el("td");
  td.colSpan = COLUMNS.length;
  const close = el("button", "close", "close");
  close.onclick = e => { e.stopPropagation(); openRank = null; renderTable(); };
  td.appendChild(close);
  const panel = el("div");
  td.appendChild(panel);
  renderPair(panel, h.model);
  tr.appendChild(td);
  tr.onclick = e => e.stopPropagation();
  return tr;
}

document.getElementById("title").textContent = `${R.query.label}: ${R.hits.length} of ${R.n_targets_tested} targets, ${R.alphabet}, k=${R.ksize}`;
document.getElementById("meta").textContent = `${R.query.name} (${R.query.length} aa). Rows are ordered by Benjamini-Hochberg corrected region tail probability (q-value); click a column to sort, click a row to open its alignments.`;
document.getElementById("qtrack").appendChild(queryTrack(R.query));
renderTable();
"""

REPORT_PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>__PAIR_CSS__ __REPORT_CSS__</style>
</head>
<body>
<h1 id="title"></h1>
<p class="meta" id="meta"></p>
<div class="legend">
  <span><i class="trk"></i> query protein, with its domains as boxes</span>
  <span><i class="seg"></i> run of 2 or more consecutive shared k-mers, on the query</span>
  <span><i class="dotm"></i> single shared k-mer</span>
</div>
<div class="qtrack" id="qtrack"></div>
<div class="scroll"><table id="hits"></table></div>
<script>
__PAIR_JS__
__REPORT_JS__
</script>
</body>
</html>
"""


def render_report(report):
    title = f"{report['query']['label']} kmerseek hits"
    return (
        REPORT_PAGE.replace("__TITLE__", html.escape(title))
        .replace("__PAIR_CSS__", PAIR_CSS)
        .replace("__REPORT_CSS__", REPORT_CSS)
        .replace("__PAIR_JS__", PAIR_JS)
        .replace("__REPORT_JS__", REPORT_JS.replace("__REPORT__", embed_json(report)))
    )


# -- CLI -----------------------------------------------------------------------------------


def _build_arg_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv", required=True, help="CSV written by `kmerseek search -o`")
    p.add_argument("--query-fasta", required=True, help="the FASTA the search was run with")
    p.add_argument("--target-fasta", required=True, help="the FASTA the target index was built from (gzip or plain)")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--query-name", help="one query's header, or its first token; default every query in the CSV")
    p.add_argument("--domains", nargs="*", default=[], metavar="TABLE", help="Pfam-style domain tables for queries and targets")
    p.add_argument("--max-hits", type=int, default=100, help="targets per query, best q-value first (default 100)")
    p.add_argument("--flank", type=int, default=0, help="residues shown either side of each run in the alignments")
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
