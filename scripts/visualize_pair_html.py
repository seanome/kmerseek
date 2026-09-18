"""Self-contained interactive HTML for one pair, and the pieces the search report reuses.

`PAIR_CSS` and `PAIR_JS` draw one figure model (see pair_model.build_model) into a
container: the dot plot with protein tracks and numbered runs, then one alignment block
per run. Hovering a single shows its k-mer; clicking a run's segment or number scrolls
to its block. `render_html(model)` wraps that in a page with the model embedded, so the
file works offline with no external libraries.
"""

import html
import json

PAIR_CSS = r"""
  :root {
    --ink: #0b0b0b; --secondary: #52514e; --surface: #fcfcfb; --edge: #b8b7b0; --box: #ffffff;
    --run: #0b0b0b; --single: #6e6d68; --domain: #e8e7e2; --domain-edge: #8d8c86; --flash: #fff3c4; --structure: #c0392b;
    --c0: #f3e3c3; --c0e: #c9a15a; --c1: #dde9f8; --c1e: #7fa6d6; --c2: #dff0e2; --c2e: #6fae7c; --c3: #ece0f3; --c3e: #a98bc4;
  }
  body { margin: 0; padding: 16px; background: var(--surface); color: var(--ink);
         font: 13px/1.45 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; }
  .pair h2 { font-size: 14px; font-weight: 600; margin: 0 0 2px; }
  .pair .sub { color: var(--secondary); margin: 0 0 10px; }
  .legend { display: flex; flex-wrap: wrap; gap: 6px 18px; color: var(--secondary); margin-bottom: 10px; max-width: 900px; }
  .legend span { display: inline-flex; align-items: center; gap: 6px; }
  .sw { width: 14px; height: 12px; display: inline-block; border: 1px solid var(--edge); box-sizing: border-box; border-radius: 2px; }
  .seg { width: 22px; height: 0; border-top: 3px solid var(--run); display: inline-block; }
  .dotm { width: 8px; height: 8px; border-radius: 50%; background: var(--single); display: inline-block; }
  .path { width: 22px; height: 0; border-top: 1.5px solid var(--structure); display: inline-block; }
  .trk { width: 26px; height: 12px; display: inline-block; position: relative; }
  .trk:before { content: ""; position: absolute; left: 0; right: 0; top: 5px; border-top: 1px solid var(--domain-edge); }
  .trk:after { content: ""; position: absolute; left: 8px; width: 10px; top: 1px; height: 9px; background: var(--domain); border: 1px solid var(--domain-edge); }
  .mid { font-family: ui-monospace, Menlo, Consolas, monospace; font-weight: 600; }
  svg text { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; font-size: 11px; fill: var(--secondary); }
  svg text.mono { font-family: ui-monospace, Menlo, Consolas, monospace; font-size: 12px; fill: var(--ink); }
  svg text.num { font-size: 13px; font-weight: 700; fill: var(--run); }
  .run-seg, .run-num { cursor: pointer; }
  .block { margin: 14px 0 0; padding: 6px 8px; border-radius: 4px; transition: background .6s; }
  .block.flash { background: var(--flash); transition: none; }
  .block h3 { font-size: 14px; font-weight: 600; margin: 0 0 4px; }
  .block .wrap { overflow-x: auto; }
  .none { color: var(--secondary); }
  #tip { position: fixed; pointer-events: none; background: #fff; border: 1px solid var(--edge); padding: 6px 8px;
         font-size: 12px; box-shadow: 0 2px 6px rgba(0,0,0,.12); display: none; white-space: nowrap; z-index: 10; }
  #tip b { font-family: ui-monospace, Menlo, Consolas, monospace; font-weight: 600; }
"""

PAIR_JS = r"""
const CLASS_FILL = ["var(--c0)", "var(--c1)", "var(--c2)", "var(--c3)"];
const CLASS_EDGE = ["var(--c0e)", "var(--c1e)", "var(--c2e)", "var(--c3e)"];
const WRAP = 60, CELL = 16, BOX = 14, ROW = 18;

const svgEl = (tag, attrs = {}, text) => {
  const e = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (text !== undefined) e.textContent = text;
  return e;
};
const el = (tag, cls, text) => { const e = document.createElement(tag); if (cls) e.className = cls; if (text !== undefined) e.textContent = text; return e; };
const range1 = (s, e) => `${s + 1}–${e}`;

function ensureTip() {
  let tip = document.getElementById("tip");
  if (!tip) { tip = el("div"); tip.id = "tip"; document.body.appendChild(tip); }
  return tip;
}

function classStyles(m) {
  const styles = {};
  if (m.classes.length && m.classes.length <= CLASS_FILL.length)
    m.classes.forEach((c, i) => styles[c.symbol] = [CLASS_FILL[i], CLASS_EDGE[i]]);
  return styles;
}

function runHeader(b) {
  const parts = [`Run ${b.number}`, `${b.query_region} × ${b.target_region}`, `${b.length} aa`, `${b.identical} identical`];
  if (b.gapped) parts[3] = `${b.identical} identical on the run's diagonal, ${b.gapped.identical} of its ${b.gapped.aligned} aligned columns after gapped alignment`;
  if (b.polar !== null && b.polar !== undefined) parts.push(`${b.polar} of ${b.length} polar`);
  if ("structure_offset" in b) parts.push(structurePhrase(b.structure_offset));
  return parts.join(" · ");
}

function blockRows(b) {
  // The gapped alignment when the model has one, else the exact run with any flank.
  if (b.gapped) return b.gapped;
  const left = b.query_start - b.window.query_start;
  return { query_row: b.query_row, target_row: b.target_row, query_enc: b.query_enc, target_enc: b.target_enc, middle: b.middle,
           query_start: b.window.query_start, target_start: b.window.target_start, run_columns: [left, left + b.length] };
}
const residuesBefore = (row, n) => [...row.slice(0, n)].filter(c => c !== "-").length;

function titleLines(m) {
  const q = m.query.label, t = m.target.label, k = m.ksize;
  const first = `${q} (query) vs ${t} (target): ${m.n_shared} shared ${k}-mers in the ${m.alphabet}`;
  const kinds = m.classes.map(c => c.label.split(" (")[0]).join("/");
  const second = m.classes.length
    ? `a shared ${k}-mer is ${k} consecutive residues with the same ${kinds} pattern in both proteins`
    : `a shared ${k}-mer is ${k} consecutive identical residues in both proteins`;
  const lines = [first, second];
  if (m.structure) {
    const st = m.structure;
    lines.push(`${st.aligner} of ${st.query_file} against ${st.target_file}: TM-score ${st.tm_score_query.toFixed(2)} (by ${q} length), RMSD ${st.rmsd.toFixed(1)} \u00c5 over ${st.aligned} aligned residues`);
  }
  return lines;
}

function structurePhrase(offset) {
  if (offset === null || offset === undefined) return "not structurally aligned";
  if (Math.abs(offset) <= 1) return "on the structural path";
  return `${Math.abs(offset)} residues off the structural path`;
}

function renderLegend(m, styles) {
  const box = el("div", "legend");
  const add = (mark, text) => { const s = el("span"); s.appendChild(mark); s.appendChild(document.createTextNode(text)); box.appendChild(s); };
  if (m.classes.length) for (const c of m.classes) {
    const sw = el("i", "sw"); const st = styles[c.symbol] || ["var(--box)", "var(--edge)"];
    sw.style.background = st[0]; sw.style.borderColor = st[1]; add(sw, c.label);
  } else add(el("i", "sw"), "residue");
  add(el("i", "seg"), `run of 2 or more consecutive shared ${m.ksize}-mers (${m.runs.length}), numbered and underlined in its alignment; click one to jump to it`);
  add(el("i", "dotm"), `single shared ${m.ksize}-mer (${m.singles.length}); hover for the k-mer`);
  if (m.query.domains.length || m.target.domains.length)
    add(el("i", "trk"), "protein, with its domains as boxes; each domain's span shaded across the plot");
  if (m.structure)
    add(el("i", "path"), `structural alignment (${m.structure.aligner}, TM-score ${m.structure.tm_score_query.toFixed(2)}): every aligned residue pair`);
  add(el("b", "mid", "G"), "identical residue, written between the rows");
  return box;
}

function pathSegments(pairs) {
  // Runs of consecutive residue pairs, so gaps break the drawn line.
  const segs = []; let cur = [];
  for (const [q, t] of pairs) {
    if (cur.length && (q !== cur[cur.length - 1][0] + 1 || t !== cur[cur.length - 1][1] + 1)) { segs.push(cur); cur = []; }
    cur.push([q, t]);
  }
  if (cur.length) segs.push(cur);
  return segs;
}

function renderStructurePath(m, svg, sx, sy) {
  if (!m.structure) return;
  for (const seg of pathSegments(m.structure.pairs)) {
    const d = seg.map(([q, t], i) => `${i ? "L" : "M"}${sx(q + 1).toFixed(1)} ${sy(t + 1).toFixed(1)}`).join("");
    svg.appendChild(svgEl("path", { d: seg.length > 1 ? d : d + "h0.1", fill: "none", stroke: "var(--structure)", "stroke-width": 1.5, "stroke-linecap": "round" }));
  }
}

function labelsCollide(side, scale) {
  const c = side.domains.map(d => [(d.start + d.end) / 2, d.name.length]);
  for (let i = 1; i < c.length; i++) if ((c[i][0] - c[i - 1][0]) * scale < (c[i][1] + c[i - 1][1]) / 2 * 6.5 + 4) return true;
  return false;
}

function renderDotPlot(m, container) {
  const side = 340, track = 36, ml = 56, mb = 40, mt = 6, mr = 8;
  const qL = m.query.length, tL = m.target.length;
  const scale = side / Math.max(qL, tL);
  const pw = qL * scale, ph = tL * scale;
  const W = ml + pw + 6 + track + 90 + mr, H = mt + track + 6 + ph + mb;
  const svg = svgEl("svg", { width: W, height: H });
  const x0 = ml, y0 = mt + track + 6;
  const sx = p => x0 + (p - 1) * scale, sy = p => y0 + (tL - p) * scale;
  const shade = attrs => svg.appendChild(svgEl("rect", { fill: "var(--domain)", opacity: 0.6, ...attrs }));
  for (const d of m.query.domains) shade({ x: sx(d.start), y: y0, width: (d.end - d.start + 1) * scale, height: ph });
  for (const d of m.target.domains) shade({ x: x0, y: sy(d.end + 1), width: pw, height: (d.end - d.start + 1) * scale });
  // Query track above the plot, target track to its right.
  const stagger = labelsCollide(m.query, scale);
  svg.appendChild(svgEl("line", { x1: sx(1), x2: sx(qL), y1: mt + track * 0.6, y2: mt + track * 0.6, stroke: "var(--domain-edge)" }));
  m.query.domains.forEach((d, i) => {
    svg.appendChild(svgEl("rect", { x: sx(d.start), y: mt + track * 0.35, width: (d.end - d.start + 1) * scale, height: track * 0.5, fill: "var(--domain)", stroke: "var(--domain-edge)" }));
    svg.appendChild(svgEl("text", { x: sx((d.start + d.end) / 2), y: mt + track * 0.28 - (stagger && i % 2 ? 12 : 0), "text-anchor": "middle" }, d.name));
  });
  const tx = x0 + pw + 6;
  svg.appendChild(svgEl("line", { x1: tx + track * 0.4, x2: tx + track * 0.4, y1: sy(1), y2: sy(tL), stroke: "var(--domain-edge)" }));
  for (const d of m.target.domains) {
    svg.appendChild(svgEl("rect", { x: tx + track * 0.15, y: sy(d.end), width: track * 0.5, height: (d.end - d.start + 1) * scale, fill: "var(--domain)", stroke: "var(--domain-edge)" }));
    svg.appendChild(svgEl("text", { x: tx + track * 0.75, y: sy((d.start + d.end) / 2) + 4 }, d.name));
  }
  // Axes with 1-based ticks.
  svg.appendChild(svgEl("line", { x1: x0, x2: x0 + pw, y1: y0 + ph, y2: y0 + ph, stroke: "var(--secondary)" }));
  svg.appendChild(svgEl("line", { x1: x0, x2: x0, y1: y0, y2: y0 + ph, stroke: "var(--secondary)" }));
  // 1, every 100 residues (200 past 600 aa, 500 past 1500), and the length; a round tick
  // within 5% of the length would overlap its label.
  const ticks = L => { const step = L <= 600 ? 100 : L <= 1500 ? 200 : 500;
    return [1, ...Array.from({ length: Math.floor((L - 1) / step) }, (_, i) => (i + 1) * step).filter(t => L - t > L * 0.05), L]; };
  for (const p of ticks(qL)) svg.appendChild(svgEl("text", { x: sx(p), y: y0 + ph + 14, "text-anchor": "middle" }, p));
  for (const p of ticks(tL)) svg.appendChild(svgEl("text", { x: x0 - 6, y: sy(p) + 4, "text-anchor": "end" }, p));
  svg.appendChild(svgEl("text", { x: x0 + pw / 2, y: H - 4, "text-anchor": "middle" }, `${m.query.label} position (aa)`));
  svg.appendChild(svgEl("text", { x: 12, y: y0 + ph / 2, "text-anchor": "middle", transform: `rotate(-90 12 ${y0 + ph / 2})` }, `${m.target.label} position (aa)`));
  renderStructurePath(m, svg, sx, sy);
  renderSingles(m, svg, sx, sy);
  renderRuns(m, svg, sx, sy, container);
  return svg;
}

function renderSingles(m, svg, sx, sy) {
  // A dot at the centre of the k residues; hover shows the k-mer.
  const tip = ensureTip(), k = m.ksize;
  for (const s of m.singles) {
    const c = svgEl("circle", { cx: sx(s.query_pos + (k + 1) / 2), cy: sy(s.target_pos + (k + 1) / 2), r: 3, fill: "var(--single)" });
    c.addEventListener("mousemove", e => {
      tip.style.display = "block";
      tip.innerHTML = `<b>${s.kmer}</b><br>${m.query.label} ${range1(s.query_pos, s.query_pos + k)}: <b>${s.query_kmer}</b><br>${m.target.label} ${range1(s.target_pos, s.target_pos + k)}: <b>${s.target_kmer}</b>`;
      tip.style.left = (e.clientX + 12) + "px"; tip.style.top = (e.clientY + 12) + "px";
    });
    c.addEventListener("mouseleave", () => tip.style.display = "none");
    svg.appendChild(c);
  }
}

function renderRuns(m, svg, sx, sy, container) {
  // A segment over the residues the run covers, numbered at its upper end.
  for (const b of m.runs) {
    const jump = () => flashBlock(container, b.number);
    const seg = svgEl("line", { x1: sx(b.query_start + 1), y1: sy(b.target_start + 1), x2: sx(b.query_end), y2: sy(b.target_end), stroke: "var(--run)", "stroke-width": 3, "stroke-linecap": "round", class: "run-seg" });
    const num = svgEl("text", { x: sx(b.query_end) + 5, y: sy(b.target_end) - 1, class: "num run-num" }, b.number);
    seg.addEventListener("click", jump); num.addEventListener("click", jump);
    svg.appendChild(seg); svg.appendChild(num);
  }
}

function flashBlock(container, number) {
  const block = container.querySelector(`[data-run="${number}"]`);
  if (!block) return;
  block.scrollIntoView({ behavior: "smooth", block: "center" });
  block.classList.add("flash");
  setTimeout(() => block.classList.remove("flash"), 900);
}

function renderBlock(m, b, styles) {
  const div = el("div", "block"); div.dataset.run = b.number;
  div.appendChild(el("h3", null, runHeader(b)));
  const wrap = el("div", "wrap"); div.appendChild(wrap);
  const rows = blockRows(b);
  for (let start = 0; start < rows.query_row.length; start += WRAP) wrap.appendChild(renderChunk(m, rows, styles, start));
  return div;
}

function renderChunk(m, rows, styles, start) {
  const sl = s => s.slice(start, start + WRAP);
  const qRow = sl(rows.query_row), tRow = sl(rows.target_row), qEnc = sl(rows.query_enc), tEnc = sl(rows.target_enc), mid = sl(rows.middle);
  const n = qRow.length;
  const x0 = Math.max(m.query.label.length, m.target.label.length) * 7 + 52;
  const svg = svgEl("svg", { width: x0 + n * CELL + 50, height: 3 * ROW + 12 });
  const row = (y, text, enc, label, full, origin) => {
    const first = origin + residuesBefore(full, start), last = origin + residuesBefore(full, start + n);
    svg.appendChild(svgEl("text", { x: x0 - 8, y: y + BOX / 2 + 4, "text-anchor": "end" }, `${label} ${first + 1}`));
    for (let i = 0; i < n; i++) {
      const st = styles[enc[i]] || ["var(--box)", "var(--edge)"];
      svg.appendChild(svgEl("rect", { x: x0 + i * CELL, y, width: BOX, height: BOX, rx: 2, fill: st[0], stroke: st[1] }));
      svg.appendChild(svgEl("text", { x: x0 + i * CELL + BOX / 2, y: y + BOX / 2 + 4, "text-anchor": "middle", class: "mono" }, text[i]));
    }
    svg.appendChild(svgEl("text", { x: x0 + n * CELL + 4, y: y + BOX / 2 + 4 }, last));
  };
  row(2, qRow, qEnc, m.query.label, rows.query_row, rows.query_start);
  for (let i = 0; i < n; i++) if (mid[i] !== " ")
    svg.appendChild(svgEl("text", { x: x0 + i * CELL + BOX / 2, y: ROW + BOX / 2 + 5, "text-anchor": "middle", class: "mono" }, mid[i]));
  row(2 * ROW + 2, tRow, tEnc, m.target.label, rows.target_row, rows.target_start);
  // The run's columns, underlined with the dot plot's run bar.
  const first = Math.max(rows.run_columns[0], start) - start, last = Math.min(rows.run_columns[1], start + n) - start;
  if (last > first) svg.appendChild(svgEl("line", { x1: x0 + first * CELL, x2: x0 + (last - 1) * CELL + BOX, y1: 3 * ROW + 7, y2: 3 * ROW + 7, stroke: "var(--run)", "stroke-width": 3 }));
  return svg;
}

function renderPair(container, m) {
  container.classList.add("pair");
  container.innerHTML = "";
  const styles = classStyles(m);
  const [first, ...rest] = titleLines(m);
  container.appendChild(el("h2", null, first));
  for (const line of rest) container.appendChild(el("p", "sub", line));
  container.appendChild(renderLegend(m, styles));
  container.appendChild(renderDotPlot(m, container));
  if (!m.runs.length) container.appendChild(el("p", "none", `No run to show: no two shared ${m.ksize}-mers are consecutive in both proteins.`));
  for (const b of m.runs) container.appendChild(renderBlock(m, b, styles));
}
"""

PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>__CSS__</style>
</head>
<body>
<div id="pair"></div>
<script>
__JS__
renderPair(document.getElementById("pair"), __MODEL__);
</script>
</body>
</html>
"""


def embed_json(value):
    """JSON safe inside a <script>: `</` cannot close the tag."""
    return json.dumps(value).replace("</", "<\\/")


def render_html(model):
    title = f"{model['query']['label']} vs {model['target']['label']} shared k-mers"
    return (
        PAGE.replace("__TITLE__", html.escape(title))
        .replace("__CSS__", PAIR_CSS)
        .replace("__JS__", PAIR_JS)
        .replace("__MODEL__", embed_json(model))
    )
