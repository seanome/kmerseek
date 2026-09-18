#!/usr/bin/env python3
"""Render one interactive HTML report per query from `kmerseek search` output.

The query is the shared axis. It is drawn once at the top as a line with its domains
(from --domains), under a histogram of how many database entries have a run over each
residue: grey for any run, black for a run with --solid-identical or more identical
residues. That histogram is the noise map: a low-complexity stretch such as the BCL-2
loop is covered by a third of unrelated proteins, so a run there is discounted at a
glance and a run in BH1 is not.

Below it, one row per protein with the numbers in the row (length, runs, longest run,
identical residues in it, shared k-mers, the ranking statistic) and every run drawn as
a bar at its query coordinates, solid when it has --solid-identical or more identical
residues and hollow otherwise; overlapping bars get a count. Database entries of one
gene (UniProt GN= and OS=) fold into one row, so a family search is not a list of
TrEMBL copies of the query. Clicking a row opens the pair view underneath it: the dot
plot with protein tracks and one alignment block per run, the same panel
visualize_pair.py draws.

Rows are ordered by `region_evalue` when the CSV has it (kmerseek >= 0.5) and otherwise
by the Benjamini-Hochberg corrected region tail probability. Sorting by identical
residues or run length instead puts composition-driven hits (p53, POU4F1) among family
members, which is why the ranking statistic is the default.

The pair view needs every shared k-mer, which the CSV does not carry, so this script runs
`kmerseek pair` once per row on sequences taken from the two FASTA files.

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
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pair_model import build_model, count_agreement, domains_for, load_domains
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


def run_bar(block):
    """What a row's bar and its tooltip need, 0-based half-open as in the model."""
    return {k: block[k] for k in ("number", "query_start", "query_end", "target_start", "target_end", "length", "identical", "polar")}


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
        "best_identical": best["identical"] if best else 0,
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
            model = build_model(pair, self.domain_rows, flank=self.args.flank)
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
            "max_runs_shown": self.args.max_runs_shown,
            "coverage": coverage(rows, query["length"], ksize, self.args.solid_identical),
            "rows": protein_rows,
        }

    def render(self, query_name, rows):
        return render_report(self.report(query_name, rows))


# -- HTML ----------------------------------------------------------------------------------

REPORT_CSS = r"""
  h1 { font-size: 16px; font-weight: 600; margin: 0 0 4px; }
  .sub { color: var(--secondary); margin: 0 0 10px; max-width: 90ch; }
  .controls { display: flex; flex-wrap: wrap; gap: 16px; align-items: center; margin: 0 0 8px; color: var(--secondary); font-size: 12px; }
  select { font: inherit; padding: 2px 4px; }
  .g { display: grid; grid-template-columns: 200px 44px 44px 64px 64px 60px 80px 300px; column-gap: 10px; align-items: center;
       padding: 4px 6px; border-bottom: 1px solid #e3e2dc; }
  .g.head { color: var(--secondary); font-size: 12px; border-bottom: 1px solid var(--edge); }
  .g.row { cursor: pointer; }
  .g.row:hover { background: #f1f0eb; }
  .g.row.open { background: #ecebe4; }
  .num { text-align: right; font-variant-numeric: tabular-nums; }
  .small { font-size: 12px; color: var(--secondary); }
  .exp { padding: 8px 6px 14px 22px; border-bottom: 1px solid #e3e2dc; background: #fafaf7; }
  .exp .shown { color: var(--secondary); font-size: 12px; margin: 6px 0 0; }
  rect.reg { fill: var(--domain); stroke: var(--domain-edge); stroke-width: 0.6; }
  rect.run { fill: var(--run); }
  rect.runlow { fill: none; stroke: var(--run); stroke-width: 1; }
  path.cov { fill: var(--single); }
  path.cov5 { fill: var(--run); }
  line.pl { stroke: var(--domain-edge); stroke-width: 1.2; }
  line.ax { stroke: var(--edge); stroke-width: 0.6; }
  text.pile { font-size: 9px; fill: var(--run); font-weight: 600; }
  .legend .runlow-sw { width: 18px; height: 8px; display: inline-block; border: 1px solid var(--run); box-sizing: border-box; }
  .legend .cov-sw { display: inline-block; width: 14px; height: 12px; position: relative; }
  .legend .cov-sw:before { content: ""; position: absolute; left: 0; width: 6px; top: 0; bottom: 0; background: var(--single); }
  .legend .cov-sw:after { content: ""; position: absolute; left: 8px; width: 6px; top: 6px; bottom: 0; background: var(--run); }
  .scroll { overflow-x: auto; }
"""

REPORT_JS = r"""
const R = __REPORT__;
const K = R.ksize, Q = R.query, QL = Q.length, TW = 300, PX = TW / QL, SOLID = R.solid_identical;
const X = p => (p - 1) * PX;   // p is 1-based
const fmtStat = v => v === 0 ? "0" : v < 1e-3 ? v.toExponential(2) : v.toFixed(3);
const bestRun = r => r.runs[0] || { length: 0, identical: 0 };
const SORTS = {
  stat: [r => [r.stat, -bestRun(r).length], `${R.stat_name} (${R.stat_note})`],
  identical: [r => [-bestRun(r).identical, -bestRun(r).length], "identical residues in the longest run"],
  longest: [r => [-bestRun(r).length, -bestRun(r).identical], "length of the longest run"],
  runs: [r => [-r.runs.length, -bestRun(r).length], "number of runs"],
  shared: [r => [-r.n_shared, -bestRun(r).length], `shared ${K}-mers`],
};
let open = new Set();

function titles() {
  const kinds = R.classes.map(c => c.label.split(" (")[0]).join("/");
  document.getElementById("title").textContent =
    `${Q.label} against ${R.n_entries} database entries (${R.n_proteins} proteins): shared ${K}-mers in the ${R.alphabet}`;
  document.getElementById("definition").textContent =
    (R.classes.length ? `A shared ${K}-mer is ${K} consecutive residues with the same ${kinds} pattern in both proteins. ` : `A shared ${K}-mer is ${K} consecutive identical residues. `) +
    `A run is 2 or more consecutive shared ${K}-mers on one diagonal. Entries of one gene fold into one row; the ` +
    `${R.rows.length} rows shown are the best by ${R.stat_name}. Click a row for its dot plot and alignments; hover a bar for its coordinates.`;
  const sel = document.getElementById("sort");
  for (const [key, [, label]] of Object.entries(SORTS)) { const o = document.createElement("option"); o.value = key; o.textContent = label; sel.appendChild(o); }
}

function reportLegend() {
  const box = document.getElementById("legend");
  const add = (mark, text) => { const s = el("span"); s.innerHTML = mark + " " + text; box.appendChild(s); };
  add(`<i class="seg"></i>`, `run of consecutive shared ${K}-mers with ${SOLID} or more identical residues`);
  add(`<i class="runlow-sw"></i>`, `run with fewer than ${SOLID} identical residues`);
  add(`<i class="cov-sw"></i>`, `database entries with any run over this query residue (grey) and with a ${SOLID}-or-more-identical run (black)`);
  add(`<i class="trk"></i>`, "query protein, with its domains as boxes");
  for (const c of R.classes) { const st = classStyles({ classes: R.classes })[c.symbol]; add(`<i class="sw" style="background:${st[0]};border-color:${st[1]}"></i>`, c.label); }
  add(`<b class="mid">G</b>`, "identical residue, written between the rows");
}

function areaPath(vals, mx, H) {
  let d = `M0 ${H}`;
  vals.forEach((v, i) => { const y = (H - v / mx * H).toFixed(1); d += `L${X(i + 1).toFixed(1)} ${y}L${X(i + 2).toFixed(1)} ${y}`; });
  return d + `L${TW} ${H}Z`;
}

function histogram() {
  const { any, solid } = R.coverage, mx = Math.max(1, ...any), H = 40;
  const svg = svgEl("svg", { width: TW, height: H + 4 });
  svg.appendChild(svgEl("path", { class: "cov", d: areaPath(any, mx, H) }));
  svg.appendChild(svgEl("path", { class: "cov5", d: areaPath(solid, mx, H) }));
  svg.appendChild(svgEl("line", { class: "ax", x1: 0, x2: TW, y1: H, y2: H }));
  return [svg, mx, any.indexOf(mx) + 1, Math.max(0, ...solid)];
}

function queryLine() {
  const svg = svgEl("svg", { width: TW, height: 34 });
  svg.appendChild(svgEl("line", { class: "pl", x1: 0, x2: TW, y1: 9, y2: 9 }));
  const stagger = labelsCollide(Q, PX);
  Q.domains.forEach((d, i) => {
    svg.appendChild(svgEl("rect", { class: "reg", x: X(d.start), y: 2, width: X(d.end + 1) - X(d.start), height: 14, rx: 2 }));
    svg.appendChild(svgEl("text", { x: (X(d.start) + X(d.end + 1)) / 2, y: 30 - (stagger && i % 2 ? 0 : 0), "text-anchor": "middle" }, d.name));
  });
  return svg;
}

function piles(runs) {
  // Runs that overlap on the query, so a count can sit over the pile.
  const sorted = [...runs].sort((a, b) => a.query_start - b.query_start), out = [];
  for (const r of sorted) {
    const last = out[out.length - 1];
    if (last && r.query_start < last.end) { last.n++; last.end = Math.max(last.end, r.query_end); }
    else out.push({ start: r.query_start, end: r.query_end, n: 1 });
  }
  return out.filter(p => p.n > 1);
}

function track(row) {
  const svg = svgEl("svg", { width: TW, height: 18 });
  for (const r of row.runs) {
    const solid = r.identical >= SOLID;
    const rect = svgEl("rect", solid
      ? { class: "run", x: X(r.query_start + 1), y: 6, width: Math.max(r.length * PX, 2), height: 8 }
      : { class: "runlow", x: X(r.query_start + 1) + 0.5, y: 6.5, width: Math.max(r.length * PX - 1, 1), height: 7 });
    rect.appendChild(svgEl("title", {}, `${Q.label} ${r.query_start + 1}–${r.query_end} × ${row.label} ${r.target_start + 1}–${r.target_end}: ${r.length} aa, ${r.identical} identical` + (r.polar === null ? "" : `, ${r.polar} polar`)));
    svg.appendChild(rect);
  }
  for (const p of piles(row.runs)) svg.appendChild(svgEl("text", { class: "pile", x: (X(p.start + 1) + X(p.end + 1)) / 2, y: 5, "text-anchor": "middle" }, p.n));
  return svg;
}

function cell(cls, text) { const d = el("div", cls); if (text !== undefined) d.textContent = text; return d; }

function headerRows(table) {
  const [hist, mx, at, mxSolid] = histogram();
  const h = el("div", "g");
  const hl = cell(null, "database entries with a run over this residue");
  hl.appendChild(cell("small", `max ${mx} of ${R.n_entries}, at residue ${at}; ${mxSolid} with a ${SOLID}-or-more-identical run`));
  h.appendChild(hl); for (let i = 0; i < 6; i++) h.appendChild(cell()); h.appendChild(hist); table.appendChild(h);
  const q = el("div", "g");
  const ql = cell(null, `${Q.label}, the query`); ql.appendChild(cell("small", `${QL} aa`));
  q.appendChild(ql); for (let i = 0; i < 6; i++) q.appendChild(cell()); q.appendChild(queryLine()); table.appendChild(q);
  const head = el("div", "g head");
  for (const [cls, text] of [[null, "target, one row per protein"], ["num", "aa"], ["num", "runs"], ["num", "longest run, aa"], ["num", "identical in it"], ["num", `shared ${K}-mers`], ["num", R.stat_name], [null, `runs drawn on the query (residue 1 to ${QL})`]])
    head.appendChild(cell(cls, text));
  table.appendChild(head);
}

function proteinRow(row) {
  const g = el("div", "g row" + (open.has(row.rank) ? " open" : ""));
  const name = cell(null, row.description || row.label);
  const extra = row.n_entries > 1 ? ` · ${row.n_entries - 1} more ${row.n_entries > 2 ? "entries" : "entry"} of this gene` : "";
  name.appendChild(cell("small", `${row.entry.split("|").pop()}${row.gene ? " · " + row.gene : ""}${extra}`));
  name.title = row.target_name + (row.other_entries.length ? "\nalso: " + row.other_entries.join(", ") : "");
  g.appendChild(name);
  for (const v of [row.length, row.runs.length, row.best_length, row.best_identical, row.n_shared, fmtStat(row.stat)]) g.appendChild(cell("num", v));
  g.appendChild(track(row));
  g.onclick = () => { open.has(row.rank) ? open.delete(row.rank) : open.add(row.rank); renderTable(); };
  return g;
}

function expansion(row) {
  const d = el("div", "exp");
  const panel = el("div");
  d.appendChild(panel);
  const shown = Math.min(row.runs.length, R.max_runs_shown);
  const m = shown < row.runs.length ? { ...row.model, runs: row.model.runs.slice(0, shown) } : row.model;
  renderPair(panel, m);
  if (shown < row.runs.length) d.appendChild(el("p", "shown", `Showing the ${shown} longest of ${row.runs.length} runs.`));
  return d;
}

function renderTable() {
  const table = document.getElementById("table");
  table.innerHTML = "";
  headerRows(table);
  const key = SORTS[document.getElementById("sort").value][0];
  const cmp = (a, b) => { const ka = key(a), kb = key(b); for (let i = 0; i < ka.length; i++) if (ka[i] !== kb[i]) return ka[i] < kb[i] ? -1 : 1; return 0; };
  for (const row of [...R.rows].sort(cmp)) {
    table.appendChild(proteinRow(row));
    if (open.has(row.rank)) table.appendChild(expansion(row));
  }
}

titles(); reportLegend();
document.getElementById("sort").onchange = renderTable;
document.getElementById("expandall").onchange = e => { open = e.target.checked ? new Set(R.rows.map(r => r.rank)) : new Set(); renderTable(); };
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
<p class="sub" id="definition"></p>
<div class="controls">
  <label>sort rows by <select id="sort"></select></label>
  <label><input type="checkbox" id="expandall"> expand every row</label>
</div>
<div class="legend" id="legend"></div>
<div class="scroll"><div id="table"></div></div>
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
    p.add_argument("--max-rows", type=int, default=100, help="protein rows per query, best ranking statistic first (default 100)")
    p.add_argument("--max-runs-shown", type=int, default=10, help="alignments per opened row, longest first (default 10)")
    p.add_argument("--solid-identical", type=int, default=5, help="identical residues from which a run's bar is drawn solid (default 5)")
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
