"""Self-contained interactive HTML for one `kmerseek pair` JSON.

Same two panels as the static figure in visualize_pair.py, drawn in the browser with no
external libraries: hovering a dot shows the k-mer and its position in both sequences,
clicking a run (or picking one from the list) moves the residue ribbon onto it, and the
flank is a live control. The pair JSON is embedded, so the file works offline.
"""

import json

# Colours match visualize_pair.py: one meaning per colour in both panels.
TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
  :root {
    --ink: #0b0b0b; --secondary: #52514e; --muted: #898781; --surface: #fcfcfb;
    --edge: #b8b7b0; --box: #ffffff; --run: #0b0b0b; --single: #898781;
    --class0: #e08a1e; --class1: #2b6cb0; --class2: #3a9d5d; --class3: #8e5bb5;
  }
  body { margin: 0; padding: 16px; background: var(--surface); color: var(--ink);
         font: 13px/1.4 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; }
  h1 { font-size: 15px; font-weight: 600; margin: 0 0 2px; }
  .sub { color: var(--secondary); margin: 0 0 14px; }
  .panel { margin-bottom: 22px; }
  .legend { display: flex; flex-wrap: wrap; gap: 14px; align-items: center;
            color: var(--secondary); margin-bottom: 6px; }
  .legend span { display: inline-flex; align-items: center; gap: 5px; }
  .swatch { width: 14px; height: 12px; display: inline-block; border: 1px solid var(--edge);
            box-sizing: border-box; }
  .dot { width: 9px; height: 9px; border-radius: 50%; display: inline-block; }
  .dash { width: 18px; height: 0; border-top: 2px dashed var(--run); display: inline-block; }
  .tick { width: 0; height: 12px; border-left: 1.5px solid var(--secondary); display: inline-block; }
  .controls { display: flex; flex-wrap: wrap; gap: 16px; align-items: center; margin-bottom: 8px; }
  .controls label { color: var(--secondary); }
  select, input[type=number] { font: inherit; padding: 2px 4px; }
  .ribbon-wrap { overflow-x: auto; }
  svg text { font-family: ui-monospace, Menlo, Consolas, monospace; }
  svg .label { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; fill: var(--secondary); }
  svg .axis { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; fill: var(--secondary); font-size: 11px; }
  .run-box { cursor: pointer; pointer-events: all; }
  .run-box:hover { stroke-width: 2; }
  #tip { position: fixed; pointer-events: none; background: #fff; border: 1px solid var(--edge);
         padding: 6px 8px; font-size: 12px; box-shadow: 0 2px 6px rgba(0,0,0,.12); display: none;
         white-space: nowrap; }
  #tip b { font-family: ui-monospace, Menlo, Consolas, monospace; font-weight: 600; }
  .note { color: var(--secondary); }
</style>
</head>
<body>
<h1 id="title"></h1>
<p class="sub" id="subtitle"></p>

<div class="panel">
  <div class="legend" id="ribbon-legend"></div>
  <div class="controls">
    <label>run <select id="run-select"></select></label>
    <label>flank <input id="flank" type="number" min="0" max="200" value="__FLANK__" style="width:4em"> residues either side</label>
  </div>
  <div class="ribbon-wrap"><svg id="ribbon"></svg></div>
  <p class="note" id="ribbon-note"></p>
</div>

<div class="panel">
  <div class="legend" id="dots-legend"></div>
  <svg id="dots"></svg>
</div>

<div id="tip"></div>

<script>
const PAIR = __PAIR_JSON__;
const CLASS_COLORS = ["var(--class0)", "var(--class1)", "var(--class2)", "var(--class3)"];
const CELL = 18, BOX = 15, ROW = 20, RIBBON_LEFT = 10, FONT = 12;

const shortLabel = (name) => {
  const parts = name.split("|");
  if (parts.length >= 8) return parts[6].split(" ")[0];
  if (parts.length >= 3) return parts[2].split(" ")[0];
  return name.split(" ")[0];
};
const K = PAIR.ksize;
const Q = PAIR.query, T = PAIR.target;
const QL = shortLabel(Q.name), TL = shortLabel(T.name);
const reduced = Q.encoded !== Q.sequence;
// A run is a region of two or more consecutive shared k-mers; lone k-mers are singles.
const runs = PAIR.regions.filter(r => r.length > K);
const onDiagonal = (k, r) =>
  r.query_start <= k.query_pos && k.query_pos <= r.query_end - K &&
  r.target_start <= k.target_pos && k.target_pos <= r.target_end - K &&
  (k.target_pos - k.query_pos) === (r.target_start - r.query_start);
const inRun = PAIR.shared_kmers.filter(k => runs.some(r => onDiagonal(k, r)));
const single = PAIR.shared_kmers.filter(k => !runs.some(r => onDiagonal(k, r)));

// {class symbol: residues seen under it}, read off the sequences.
const classes = {};
if (reduced) for (const s of [Q, T]) for (let i = 0; i < s.sequence.length; i++) {
  (classes[s.encoded[i]] ??= new Set()).add(s.sequence[i]);
}
const classSymbols = Object.keys(classes).sort();
const classFill = {};
if (classSymbols.length <= CLASS_COLORS.length) classSymbols.forEach((c, i) => classFill[c] = CLASS_COLORS[i]);

const state = { run: 0, flank: __FLANK__ };

const svgEl = (tag, attrs = {}, text) => {
  const el = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v);
  if (text !== undefined) el.textContent = text;
  return el;
};
const agree = (a, b) => { let n = 0; for (let i = 0; i < a.length; i++) if (a[i] === b[i]) n++; return n; };
const range1 = (start, end) => `${start + 1}–${end}`;

function titles() {
  document.getElementById("title").textContent =
    `${QL} (query) vs ${TL} (target), ${PAIR.moltype} ${K}-mers`;
  const n = PAIR.shared_kmers.length;
  let sub;
  if (!runs.length) sub = `${n} shared, none consecutive in both sequences`;
  else {
    const r = runs[0];
    const ident = agree(Q.sequence.slice(r.query_start, r.query_end), T.sequence.slice(r.target_start, r.target_end));
    const same = agree(Q.encoded.slice(r.query_start, r.query_end), T.encoded.slice(r.target_start, r.target_end));
    sub = `${n} shared: ${inRun.length} in runs of consecutive k-mers, ${single.length} singles; ` +
      `longest run ${r.length} residues, ${ident}/${r.length} identical, ${same}/${r.length} same class`;
  }
  document.getElementById("subtitle").textContent = sub;
}

function ribbonLegend() {
  const el = document.getElementById("ribbon-legend");
  el.innerHTML = "";
  const add = (swatch, text) => { const s = document.createElement("span"); s.innerHTML = swatch + " " + text; el.appendChild(s); };
  add(`<i class="swatch" style="background:var(--box)"></i>`, "residue");
  for (const c of classSymbols) {
    const fill = classFill[c] ?? "var(--box)";
    add(`<i class="swatch" style="background:${fill};border-color:${classFill[c] ? fill : "var(--edge)"}"></i>`,
        `class ${c}: ${[...classes[c]].sort().join(" ")}`);
  }
  add(`<i class="tick"></i>`, reduced ? "same class in both" : "same residue in both");
  add(`<i class="dash"></i>`, "selected run of consecutive shared k-mers");
}

function runSelect() {
  const sel = document.getElementById("run-select");
  sel.innerHTML = "";
  runs.forEach((r, i) => {
    const opt = document.createElement("option");
    opt.value = i;
    opt.textContent = `${r.length - K + 1} k-mers: ${QL} ${range1(r.query_start, r.query_end)}, ${TL} ${range1(r.target_start, r.target_end)}`;
    sel.appendChild(opt);
  });
  sel.disabled = !runs.length;
  sel.onchange = () => { state.run = +sel.value; drawRibbon(); drawDots(); };
  document.getElementById("flank").oninput = (e) => { state.flank = Math.max(0, +e.target.value || 0); drawRibbon(); };
}

function windowFor(r) {
  const left = Math.min(state.flank, r.query_start, r.target_start);
  const right = Math.min(state.flank, Q.sequence.length - r.query_end, T.sequence.length - r.target_end);
  return { qs: r.query_start - left, qe: r.query_end + right, ts: r.target_start - left, te: r.target_end + right };
}

function drawRibbon() {
  const svg = document.getElementById("ribbon");
  svg.innerHTML = "";
  const note = document.getElementById("ribbon-note");
  if (!runs.length) {
    svg.setAttribute("width", 0); svg.setAttribute("height", 0);
    note.textContent = `No run to show: no two shared ${K}-mers are consecutive in both sequences.`;
    return;
  }
  const r = runs[state.run];
  const w = windowFor(r);
  const rows = reduced
    ? [[T.sequence.slice(w.ts, w.te), false, `${TL} ${range1(w.ts, w.te)}`],
       [T.encoded.slice(w.ts, w.te), true, PAIR.moltype],
       [Q.encoded.slice(w.qs, w.qe), true, PAIR.moltype],
       [Q.sequence.slice(w.qs, w.qe), false, `${QL} ${range1(w.qs, w.qe)}`]]
    : [[T.sequence.slice(w.ts, w.te), false, `${TL} ${range1(w.ts, w.te)}`],
       [Q.sequence.slice(w.qs, w.qe), false, `${QL} ${range1(w.qs, w.qe)}`]];
  const n = w.qe - w.qs;
  const labelW = Math.max(...rows.map(x => x[2].length)) * 7 + 12;
  const x0 = RIBBON_LEFT + labelW;
  const height = rows.length * ROW + 30;
  svg.setAttribute("width", x0 + n * CELL + 10);
  svg.setAttribute("height", height);

  rows.forEach(([text, encoded, label], ri) => {
    const y = 6 + ri * ROW;
    svg.appendChild(svgEl("text", { x: x0 - 8, y: y + BOX / 2 + 4, "text-anchor": "end", class: "label", "font-size": FONT }, label));
    for (let i = 0; i < text.length; i++) {
      const fill = encoded ? (classFill[text[i]] ?? "var(--box)") : "var(--box)";
      const colored = fill !== "var(--box)";
      svg.appendChild(svgEl("rect", { x: x0 + i * CELL, y, width: BOX, height: BOX, fill,
        stroke: colored ? fill : "var(--edge)", "stroke-width": 0.8 }));
      svg.appendChild(svgEl("text", { x: x0 + i * CELL + BOX / 2, y: y + BOX / 2 + 4, "text-anchor": "middle",
        "font-size": FONT, fill: colored ? "#fff" : "var(--ink)" }, text[i]));
    }
  });
  // Ticks between the two middle rows wherever they agree.
  const top = rows[rows.length / 2 - 1][0], bottom = rows[rows.length / 2][0];
  const yTick = 6 + (rows.length / 2 - 1) * ROW + BOX;
  for (let i = 0; i < top.length; i++) if (top[i] === bottom[i]) {
    svg.appendChild(svgEl("line", { x1: x0 + i * CELL + BOX / 2, x2: x0 + i * CELL + BOX / 2, y1: yTick, y2: yTick + (ROW - BOX),
      stroke: "var(--secondary)", "stroke-width": 1.2 }));
  }
  // The selected run, outlined.
  const start = r.query_start - w.qs;
  svg.appendChild(svgEl("rect", { x: x0 + start * CELL - 2, y: 3, width: r.length * CELL - (CELL - BOX) + 4,
    height: rows.length * ROW - (ROW - BOX) + 6, fill: "none", stroke: "var(--run)", "stroke-width": 1.4, "stroke-dasharray": "5 3" }));
  const ident = agree(Q.sequence.slice(r.query_start, r.query_end), T.sequence.slice(r.target_start, r.target_end));
  const same = agree(Q.encoded.slice(r.query_start, r.query_end), T.encoded.slice(r.target_start, r.target_end));
  svg.appendChild(svgEl("text", { x: x0 + start * CELL - 2, y: rows.length * ROW + 20, class: "label", "font-size": FONT },
    `${r.length - K + 1} consecutive shared ${K}-mers, ${r.length} residues: ${ident}/${r.length} identical, ` +
    `${same}/${r.length} same class` + (reduced ? "" : "")));
  note.textContent = "";
}

function dotsLegend() {
  const el = document.getElementById("dots-legend");
  el.innerHTML = "";
  const add = (swatch, text) => { const s = document.createElement("span"); s.innerHTML = swatch + " " + text; el.appendChild(s); };
  add(`<i class="dot" style="background:var(--single)"></i>`, `single shared ${K}-mer (${single.length})`);
  add(`<i class="dot" style="background:var(--run)"></i>`, `shared ${K}-mer in a run (${inRun.length})`);
  add(`<i class="dash"></i>`, `run of consecutive shared k-mers (${runs.length}); click one to show it above`);
  add(``, `hover a dot for the k-mer`);
}

function drawDots() {
  const svg = document.getElementById("dots");
  svg.innerHTML = "";
  const side = 360, ml = 58, mb = 44, mt = 8, mr = 12;
  const qLen = Q.sequence.length, tLen = T.sequence.length;
  const scale = side / Math.max(qLen, tLen);
  const W = ml + qLen * scale + mr, H = mt + tLen * scale + mb;
  svg.setAttribute("width", W); svg.setAttribute("height", H);
  const sx = (p) => ml + p * scale, sy = (p) => mt + (tLen - p) * scale;
  // Axes with residue ticks.
  svg.appendChild(svgEl("line", { x1: ml, x2: ml + qLen * scale, y1: sy(0), y2: sy(0), stroke: "var(--secondary)" }));
  svg.appendChild(svgEl("line", { x1: ml, x2: ml, y1: sy(0), y2: sy(tLen), stroke: "var(--secondary)" }));
  const step = Math.max(qLen, tLen) > 400 ? 100 : 50;
  for (let p = 0; p <= qLen; p += step) {
    svg.appendChild(svgEl("line", { x1: sx(p), x2: sx(p), y1: sy(0), y2: sy(0) + 4, stroke: "var(--secondary)" }));
    svg.appendChild(svgEl("text", { x: sx(p), y: sy(0) + 16, "text-anchor": "middle", class: "axis" }, p));
  }
  for (let p = 0; p <= tLen; p += step) {
    svg.appendChild(svgEl("line", { x1: ml - 4, x2: ml, y1: sy(p), y2: sy(p), stroke: "var(--secondary)" }));
    svg.appendChild(svgEl("text", { x: ml - 7, y: sy(p) + 4, "text-anchor": "end", class: "axis" }, p));
  }
  svg.appendChild(svgEl("text", { x: ml + qLen * scale / 2, y: H - 6, "text-anchor": "middle", class: "axis" },
    `${QL} position (query, ${qLen} aa)`));
  const yl = svgEl("text", { x: 14, y: mt + tLen * scale / 2, "text-anchor": "middle", class: "axis",
    transform: `rotate(-90 14 ${mt + tLen * scale / 2})` }, `${TL} position (target, ${tLen} aa)`);
  svg.appendChild(yl);

  const tip = document.getElementById("tip");
  const show = (e, k) => {
    tip.style.display = "block";
    tip.innerHTML = `<b>${k.kmer}</b><br>${QL} ${range1(k.query_pos, k.query_pos + K)}: <b>${k.query_kmer}</b>` +
      `<br>${TL} ${range1(k.target_pos, k.target_pos + K)}: <b>${k.target_kmer}</b>`;
    tip.style.left = (e.clientX + 12) + "px"; tip.style.top = (e.clientY + 12) + "px";
  };
  const hide = () => tip.style.display = "none";
  // A k-mer covers k residues; its dot sits half a residue past its start so position 0 is inside the axes.
  for (const [list, color] of [[single, "var(--single)"], [inRun, "var(--run)"]]) for (const k of list) {
    const c = svgEl("circle", { cx: sx(k.query_pos + 0.5), cy: sy(k.target_pos + 0.5), r: 3.2, fill: color });
    c.addEventListener("mousemove", (e) => show(e, k));
    c.addEventListener("mouseleave", hide);
    svg.appendChild(c);
  }
  runs.forEach((r, i) => {
    const box = svgEl("rect", { x: sx(r.query_start), y: sy(r.target_end), width: r.length * scale, height: r.length * scale,
      fill: "none", stroke: "var(--run)", "stroke-width": i === state.run ? 2 : 1, "stroke-dasharray": "5 3", class: "run-box" });
    box.addEventListener("click", () => { state.run = i; document.getElementById("run-select").value = i; drawRibbon(); drawDots(); });
    svg.appendChild(box);
  });
}

titles(); ribbonLegend(); runSelect(); drawRibbon(); dotsLegend(); drawDots();
</script>
</body>
</html>
"""


def render_html(pair, title, flank=10):
    """The finished page as a string. `</` inside the embedded JSON is escaped so a
    sequence name can never close the script tag."""
    pair_json = json.dumps(pair).replace("</", "<\\/")
    return (
        TEMPLATE.replace("__TITLE__", title)
        .replace("__PAIR_JSON__", pair_json)
        .replace("__FLANK__", str(flank))
    )
