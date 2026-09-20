"""The figure model for one query against one target: what both the matplotlib figure and
the HTML page draw from, so the two cannot drift.

Built from the JSON that `kmerseek pair` writes plus optional domain tables. Positions inside the model are 0-based
half-open, as in the JSON; text meant for a reader is 1-based inclusive.
"""

import os
import statistics
from collections import defaultdict

import polars as pl

from gapped_alignment import global_align
from visualize_hits import short_label

# What a biologist calls the classes of the hydrophobic/polar alphabets.
CLASS_NAMES = {"h": "hydrophobic", "p": "polar", "c": "cysteine"}
NO_REGION = "no region"

# -- pair JSON -----------------------------------------------------------------------


def is_reduced(pair):
    return pair["query"]["encoded"] != pair["query"]["sequence"]


def class_residues(pair):
    """{class symbol: sorted residues that map to it}, empty when the alphabet does not
    reduce. The alphabet's whole table when the JSON carries one (`classes`), so every
    residue is named even if these two sequences lack some; otherwise read off the two
    sequences."""
    if not is_reduced(pair):
        return {}
    if "classes" in pair:
        return {cls: "".join(sorted(res)) for cls, res in sorted(pair["classes"].items())}
    seen = defaultdict(set)
    for side in ("query", "target"):
        for residue, cls in zip(pair[side]["sequence"], pair[side]["encoded"]):
            seen[cls].add(residue)
    return {cls: "".join(sorted(res)) for cls, res in sorted(seen.items())}


def runs(pair):
    """Matched regions made of at least two consecutive shared k-mers. `kmerseek pair` also
    reports every lone shared k-mer as a region exactly k residues long; those are drawn as
    singles, not as alignments."""
    return [r for r in pair["regions"] if r["length"] > pair["ksize"]]


def in_region(kmer, region, ksize):
    """A shared k-mer sits on a region's diagonal when its start is inside the region in
    both sequences and its diagonal offset matches."""
    q_in = region["query_start"] <= kmer["query_pos"] <= region["query_end"] - ksize
    t_in = region["target_start"] <= kmer["target_pos"] <= region["target_end"] - ksize
    same_diagonal = (kmer["target_pos"] - kmer["query_pos"]) == (
        region["target_start"] - region["query_start"]
    )
    return q_in and t_in and same_diagonal


def singles(pair):
    """Shared k-mers on no run's diagonal."""
    ksize = pair["ksize"]
    return [k for k in pair["shared_kmers"] if not any(in_region(k, r, ksize) for r in runs(pair))]


def count_agreement(a, b):
    return sum(x == y for x, y in zip(a, b))


# -- domain tables ----------------------------------------------------------------------

PROTEIN_COLUMNS = ("accession", "protein", "name")
START_COLUMNS = ("domain_start", "start")
END_COLUMNS = ("domain_end", "end")
NAME_COLUMNS = ("name", "pfam_name", "pfam_id")


def _pick(columns, candidates, what):
    for c in candidates:
        if c in columns:
            return c
    raise ValueError(f"no {what} column among {list(columns)}; expected one of {candidates}")


def _read_table(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".parquet", ".pq"):
        return pl.read_parquet(path)
    return pl.read_csv(path, separator="\t" if ext in (".tsv", ".txt") else ",")


def read_domain_table(path):
    """One table of domains as rows of {protein, start, end, name}, start/end 1-based
    inclusive. Accepts the analysis repo's `*_pfam_domains.parquet` as is."""
    df = _read_table(path)
    protein = _pick(df.columns, PROTEIN_COLUMNS, "protein")
    start, end = _pick(df.columns, START_COLUMNS, "start"), _pick(df.columns, END_COLUMNS, "end")
    # `name` may be the protein column; the domain name then comes from the next candidate.
    name = _pick([c for c in df.columns if c != protein], NAME_COLUMNS, "domain name")
    rows = df.select(protein, start, end, name).drop_nulls().iter_rows(named=True)
    return [{"protein": str(r[protein]), "start": int(r[start]), "end": int(r[end]), "name": str(r[name])} for r in rows]


def load_domains(paths):
    return [row for p in paths for row in read_domain_table(p)]


def header_keys(header):
    """The strings a domain table may key this protein by: full header, first token, and
    the accession and entry name in `sp|P10415|BCL2_HUMAN ...`."""
    first = header.split()[0] if header.split() else header
    keys = {header, first}
    parts = first.split("|")
    if len(parts) >= 3:
        keys.update((parts[1], parts[2]))
    return keys


def _is_accession(name):
    """A bare Pfam accession such as PF00452, which a table without names falls back to."""
    return name.startswith("PF") and name[2:].split(".")[0].isdigit()


def domains_for(domain_rows, header):
    """This protein's domains, sorted by start. The same span listed by two tables (one
    with names, one with accessions) is kept once, under the name."""
    keys = header_keys(header)
    by_span = {}
    for d in domain_rows:
        if d["protein"] not in keys:
            continue
        span = (d["start"], d["end"])
        if span not in by_span or (_is_accession(by_span[span]["name"]) and not _is_accession(d["name"])):
            by_span[span] = {"name": d["name"], "start": d["start"], "end": d["end"]}
    return sorted(by_span.values(), key=lambda d: d["start"])


def region_name(domains, start, end):
    """The domain overlapping most of the 1-based inclusive span start..end, or NO_REGION."""
    best, best_overlap = NO_REGION, 0
    for d in domains:
        overlap = min(end, d["end"]) - max(start, d["start"]) + 1
        if overlap > best_overlap:
            best, best_overlap = d["name"], overlap
    return best


# -- the figure model ---------------------------------------------------------------------


def class_label(cls, residues):
    return f"{CLASS_NAMES.get(cls, f'class {cls}')} ({' '.join(residues)})"


def alphabet_label(moltype, classes):
    """`2-letter hydrophobic/polar alphabet (Lehninger)` for hp_lehninger2, else the
    moltype and its class count."""
    n = len(classes)
    if moltype.startswith("hp_"):
        source = moltype[3:].rstrip("0123456789").replace("_", " ").title()
        kind = "/".join(CLASS_NAMES.get(c, c) for c in classes)
        return f"{n}-letter {kind} alphabet ({source})"
    return f"{moltype} ({n} classes)" if n else f"{moltype} (no reduction)"


def middle_line(query_row, target_row, query_enc, target_enc, show_class):
    """The line between the rows: the letter where the residues are identical, `:` where
    only the class agrees (when show_class), a space otherwise."""
    out = []
    for q, t, qc, tc in zip(query_row, target_row, query_enc, target_enc):
        out.append(q if q == t else (":" if show_class and qc == tc else " "))
    return "".join(out)


def gapped_block(pair, run, gap_flank):
    """The end-to-end gapped alignment of the run and `gap_flank` residues either side, so the
    run is always traversed. Rows carry `-` for gaps; the class rows carry `-` there too, so a
    gap column is drawn as a plain box."""
    q, t = pair["query"], pair["target"]
    qs, qe, ts, te = run["query_start"], run["query_end"], run["target_start"], run["target_end"]
    q0, t0 = max(0, qs - gap_flank), max(0, ts - gap_flank)
    hit = global_align(q["sequence"][q0 : qe + gap_flank], t["sequence"][t0 : te + gap_flank])
    if hit is None:
        return None
    q_enc = _gapped_classes(hit["a_row"], q["encoded"], q0 + hit["a_start"])
    t_enc = _gapped_classes(hit["b_row"], t["encoded"], t0 + hit["b_start"])
    columns = run_columns(hit["a_row"], q0 + hit["a_start"], qs, qe)
    # Identities are counted over the run's own columns, not the flanks, so a chance match
    # in a flank cannot lift a run; the flanks are there for the gaps that shift the run.
    in_run = [(a, b) for a, b in zip(hit["a_row"][slice(*columns)], hit["b_row"][slice(*columns)])]
    return {
        "run_columns": columns,
        "query_start": q0 + hit["a_start"],
        "query_end": q0 + hit["a_end"],
        "target_start": t0 + hit["b_start"],
        "target_end": t0 + hit["b_end"],
        "query_row": hit["a_row"],
        "target_row": hit["b_row"],
        "query_enc": q_enc,
        "target_enc": t_enc,
        "middle": middle_line(hit["a_row"], hit["b_row"], q_enc, t_enc, show_class=True),
        "identical": sum(a == b and a != "-" for a, b in in_run),
        "aligned": sum(a != "-" and b != "-" for a, b in in_run),
        "window_identical": hit["identical"],
        "window_aligned": hit["aligned"],
        "columns": hit["columns"],
        "score": hit["score"],
    }


def run_columns(query_row, query_start, run_start, run_end):
    """[first, last) alignment columns whose query residue lies inside the run."""
    first, last, pos = None, None, query_start
    for col, ch in enumerate(query_row):
        if ch == "-":
            continue
        if run_start <= pos < run_end:
            first = col if first is None else first
            last = col + 1
        pos += 1
    return [first, last] if first is not None else [0, 0]


def _gapped_classes(row, encoded, start):
    """The class symbol under each column of an aligned row, `-` at a gap."""
    out, pos = [], start
    for ch in row:
        if ch == "-":
            out.append("-")
        else:
            out.append(encoded[pos])
            pos += 1
    return "".join(out)


def structural_offset(run, pairs):
    """How far the run's diagonal sits from the structural alignment over the run's query
    residues: the median of (structural target partner - run's target partner), or None when
    the structure aligns none of those residues."""
    partner = {q: t for q, t, _ in pairs}
    diagonal = run["target_start"] - run["query_start"]
    offsets = [partner[i] - (i + diagonal) for i in range(run["query_start"], run["query_end"]) if i in partner]
    return int(statistics.median(offsets)) if offsets else None


def run_block(pair, index, run, domains, flank, gap_flank=None, pairs=None):
    """Everything one alignment block needs. `structure_offset` is present only when a
    structural alignment was given."""
    q, t, ksize = pair["query"], pair["target"], pair["ksize"]
    qs, qe, ts, te = run["query_start"], run["query_end"], run["target_start"], run["target_end"]
    left = min(flank, qs, ts)
    right = min(flank, len(q["sequence"]) - qe, len(t["sequence"]) - te)
    ws, we, wts, wte = qs - left, qe + right, ts - left, te + right
    q_row, t_row, q_enc, t_enc = q["sequence"][ws:we], t["sequence"][wts:wte], q["encoded"][ws:we], t["encoded"][wts:wte]
    block = {
        "number": index + 1,
        "query_start": qs,
        "query_end": qe,
        "target_start": ts,
        "target_end": te,
        "length": run["length"],
        "n_kmers": run["length"] - ksize + 1,
        "query_region": region_name(domains["query"], qs + 1, qe),
        "target_region": region_name(domains["target"], ts + 1, te),
        "identical": count_agreement(q["sequence"][qs:qe], t["sequence"][ts:te]),
        "polar": q["encoded"][qs:qe].count("p") if "p" in class_residues(pair) else None,
        "window": {"query_start": ws, "query_end": we, "target_start": wts, "target_end": wte},
        "query_row": q_row,
        "target_row": t_row,
        "query_enc": q_enc,
        "target_enc": t_enc,
        "middle": middle_line(q_row, t_row, q_enc, t_enc, show_class=flank > 0),
        "gapped": gapped_block(pair, run, gap_flank) if gap_flank is not None else None,
    }
    if pairs is not None:
        block["structure_offset"] = structural_offset(run, pairs)
    return block


def _side(pair, side, domains):
    return {
        "label": short_label(pair[side]["name"]),
        "name": pair[side]["name"],
        "length": len(pair[side]["sequence"]),
        "domains": domains[side],
    }


def build_model(pair, domain_rows=(), flank=0, gap_flank=None, structure=None):
    """The one description both renderers draw from. `gap_flank` adds an end-to-end gapped
    alignment of each run and that many residues either side; `structure` is the dict
    `structure_alignment.align_pair` returns, or None."""
    classes = class_residues(pair)
    domains = {side: domains_for(domain_rows, pair[side]["name"]) for side in ("query", "target")}
    pairs = structure["pairs"] if structure else None
    return {
        "ksize": pair["ksize"],
        "moltype": pair["moltype"],
        "alphabet": alphabet_label(pair["moltype"], classes),
        "classes": [{"symbol": c, "residues": r, "label": class_label(c, r)} for c, r in classes.items()],
        "query": _side(pair, "query", domains),
        "target": _side(pair, "target", domains),
        "n_shared": len(pair["shared_kmers"]),
        "runs": [run_block(pair, i, r, domains, flank, gap_flank, pairs) for i, r in enumerate(runs(pair))],
        "singles": [
            {k: s[k] for k in ("query_pos", "target_pos", "kmer", "query_kmer", "target_kmer")} for s in singles(pair)
        ],
        "flank": flank,
        "gap_flank": gap_flank,
        "structure": structure,
    }


def title_lines(model):
    q, t, k = model["query"]["label"], model["target"]["label"], model["ksize"]
    first = f"{q} (query) vs {t} (target): {model['n_shared']} shared {k}-mers in the {model['alphabet']}"
    if model["classes"]:
        kinds = "/".join(CLASS_NAMES.get(c["symbol"], c["symbol"]) for c in model["classes"])
        second = f"a shared {k}-mer is {k} consecutive residues with the same {kinds} pattern in both proteins"
    else:
        second = f"a shared {k}-mer is {k} consecutive identical residues in both proteins"
    lines = [first, second]
    if model.get("structure"):
        lines.append(structure_line(model))
    return lines


def structure_line(model):
    st = model["structure"]
    return (
        f"{st['aligner']} of {st['query_file']} against {st['target_file']}: TM-score {st['tm_score_query']:.2f} "
        f"(by {model['query']['label']} length), RMSD {st['rmsd']:.1f} \u00c5 over {st['aligned']} aligned residues"
    )


def structure_phrase(offset):
    if offset is None:
        return "not structurally aligned"
    if abs(offset) <= 1:
        return "on the structural path"
    return f"{abs(offset)} residues off the structural path"


def run_header(block):
    parts = [
        f"Run {block['number']}",
        f"{block['query_region']} × {block['target_region']}",
        f"{block['length']} aa",
        f"{block['identical']} identical",
    ]
    if block.get("gapped"):
        g = block["gapped"]
        parts[-1] = f"{block['identical']} identical on the run's diagonal, {g['identical']} of its {g['aligned']} aligned columns after gapped alignment"
    if block["polar"] is not None:
        parts.append(f"{block['polar']} of {block['length']} polar")
    if "structure_offset" in block:
        parts.append(structure_phrase(block["structure_offset"]))
    return " · ".join(parts)
