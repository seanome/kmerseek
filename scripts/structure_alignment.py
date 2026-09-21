"""Superposition of a pair's structures with USalign or TM-align, as residue pairs for the dot plot.

Given a directory of structures (AlphaFold `AF-{accession}-F1-model_v*.cif` or `.pdb`, or
`{accession}.pdb`, or `{first token}.pdb`) and an aligner binary, `align_pair` runs the
aligner on the two proteins' files and parses its report: TM-scores normalised by each
length, RMSD, the aligned length, and the three alignment lines, from which every aligned
residue pair is read (0-based positions in each sequence, and whether the pair is within
5 A, which USalign marks with `:`). Drawn across the dot plot, those pairs are what the
runs get judged against: a run on them is a real correspondence, a run off them is not.
"""

import glob
import os
import re
import shutil
import subprocess

from pair_model import header_keys

ALIGNER_NAMES = ("USalign", "TMalign")


def find_aligner(explicit=None):
    """The aligner binary: --aligner, else USalign or TMalign on PATH. None when absent."""
    if explicit:
        return explicit if os.path.exists(explicit) else shutil.which(explicit)
    for name in ALIGNER_NAMES:
        found = shutil.which(name)
        if found:
            return found
    return None


def find_structure(directory, header):
    """The structure file for a protein, by accession or first header token, or None."""
    for key in sorted(header_keys(header), key=len):
        for pattern in (f"AF-{key}-F1-model_v*.cif", f"AF-{key}-F1-model_v*.pdb", f"{key}.cif", f"{key}.pdb"):
            hits = sorted(glob.glob(os.path.join(directory, pattern)))
            if hits:
                return hits[-1]  # the newest AlphaFold model version sorts last
    return None


def parse_report(text):
    """TM-scores, RMSD, aligned length and residue pairs from USalign/TM-align output."""
    tm = re.findall(r"TM-score=\s*([0-9.]+) \(normalized by length of Structure_(\d)", text)
    scores = {f"tm_score_{'query' if which == '1' else 'target'}": float(v) for v, which in tm}
    m = re.search(r"Aligned length=\s*(\d+), RMSD=\s*([0-9.]+), Seq_ID=n_identical/n_aligned=\s*([0-9.]+)", text)
    lines = [l for l in text.splitlines()]
    marks_at = next(i for i, l in enumerate(lines) if l.startswith('(":" denotes'))
    rows = [l for l in lines[marks_at + 1 :] if l and not l.startswith("#")][:3]
    return {
        **scores,
        "aligned": int(m.group(1)),
        "rmsd": float(m.group(2)),
        "seq_id": float(m.group(3)),
        "pairs": residue_pairs(rows[0], rows[1], rows[2]),
    }


def residue_pairs(row_a, marks, row_b):
    """[(a_pos, b_pos, close)] for every aligned column, 0-based; close = within 5 A."""
    pairs, i, j = [], 0, 0
    for ca, mark, cb in zip(row_a, marks.ljust(len(row_a)), row_b):
        if ca != "-" and cb != "-":
            pairs.append((i, j, mark == ":"))
        i += ca != "-"
        j += cb != "-"
    return pairs


def run_aligner(aligner, query_file, target_file):
    out = subprocess.run([aligner, query_file, target_file, "-outfmt", "0"], capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(f"{os.path.basename(aligner)} failed on {query_file} vs {target_file}:\n{out.stderr.strip()}")
    return out.stdout


def align_pair(aligner, directory, query_header, target_header):
    """The superposition for the pair, or None when either structure is missing."""
    query_file, target_file = find_structure(directory, query_header), find_structure(directory, target_header)
    if not (aligner and query_file and target_file):
        return None
    report = parse_report(run_aligner(aligner, query_file, target_file))
    report["aligner"] = os.path.basename(aligner)
    report["query_file"] = os.path.basename(query_file)
    report["target_file"] = os.path.basename(target_file)
    return report
