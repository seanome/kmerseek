#!/usr/bin/env python3
"""Draw docs/images/massspec18_classes.png: the 20 residues and the class each falls in,
under protein20, massspec18 and uniprot18, with a human peptide read both ways.

The peptide is residues 19-32 of human serum albumin (UniProt P02768), RGVFRRDAHKSEVA. A
mass spectrometer reading K28 as Q gives RGVFRRDAHQSEVA; under massspec18
both strings encode the same, under the other two they do not.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).parent))
import pubfig as pf  # noqa: E402

RESIDUES = "ACDEFGHIKLMNPQRSTVWY"
ALPHABETS = {
    "protein20": list(RESIDUES),
    "massspec18": ["A", "C", "D", "E", "F", "G", "H", "IL", "QK", "M", "N", "P", "R", "S", "T",
                   "V", "W", "Y"],
    "uniprot18": ["A", "R", "N", "D", "C", "Q", "EP", "G", "HL", "I", "K", "M", "F", "S", "T",
                  "W", "Y", "V"],
}
PAIR_COLOURS = [pf.OKABE_ITO["vermillion"], pf.OKABE_ITO["blue"]]
OWN = "#dddddd"
PEPTIDE = "RGVFRRDAHKSEVA"
READ_AS = "RGVFRRDAHQSEVA"


def encode(seq, clusters):
    symbol = {aa: c[0].lower() for c in clusters for aa in c}
    return "".join(symbol[aa] for aa in seq)


def main():
    pf.use_style()
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(pf.TWO_COLUMN_MM * pf.MM, 70 * pf.MM),
                                 layout="constrained", height_ratios=[3, 2.4])
    ax.add_patch(Rectangle((0, 3.6), 0.4, 0.3, color=PAIR_COLOURS[0], clip_on=False))
    ax.add_patch(Rectangle((0.45, 3.6), 0.4, 0.3, color=PAIR_COLOURS[1], clip_on=False))
    ax.text(0.95, 3.75, "same colour in a row = one class", va="center")
    ax.add_patch(Rectangle((7, 3.6), 0.4, 0.3, color=OWN, clip_on=False))
    ax.text(7.5, 3.75, "residue in a class of its own", va="center")
    for row, (name, clusters) in enumerate(ALPHABETS.items()):
        y = 2 - row
        cls = {aa: c for c in clusters for aa in c}
        pairs = [c for c in clusters if len(c) > 1]
        colour = {c: PAIR_COLOURS[i] for i, c in enumerate(pairs)}
        for x, aa in enumerate(RESIDUES):
            merged = len(cls[aa]) > 1
            ax.add_patch(Rectangle((x + 0.05, y + 0.1), 0.9, 0.8,
                                   color=colour[cls[aa]] if merged else OWN))
            ax.text(x + 0.5, y + 0.5, aa, ha="center", va="center",
                    color="white" if merged else "black", family="Courier New", fontsize=7)
        ax.text(-0.3, y + 0.5, f"{name} ({len(clusters)} classes)", ha="right", va="center")
    ax.set_xlim(-6.5, 20)
    ax.set_ylim(0, 4)
    ax.axis("off")
    ax.set_title("massspec18 merges only I with L and Q with K, the pairs a mass spectrometer "
                 "confuses", loc="left")

    lines = [("albumin 19-32", PEPTIDE), ("read as Q", READ_AS)]
    bx.axis("off")
    bx.set_xlim(-6.5, 20)
    bx.set_ylim(0, 6.2)
    y = 5.6
    for name in ["protein20", "massspec18", "uniprot18"]:
        a, b = (encode(s, ALPHABETS[name]) for _, s in lines)
        match = "".join("|" if p == q else " " for p, q in zip(a, b))
        n_same = match.count("|")
        bx.text(-0.3, y, f"{name}", ha="right", va="top")
        bx.text(0, y, f"{a}\n{match}\n{b}", va="top", family="Courier New", fontsize=6.5)
        verdict = "same string" if a == b else f"{n_same} of {len(a)} the same: one mismatch"
        bx.text(9, y - 0.55, verdict, va="center")
        y -= 1.9
    bx.text(-6.5, 6.2, f"Human serum albumin (P02768) residues 19-32, {PEPTIDE}, with K28 read "
            f"as Q, encoded by each alphabet", va="bottom")
    out = Path(__file__).resolve().parents[1] / "docs" / "images" / "massspec18_classes"
    pf.save(fig, out, formats=("png",))
    print(out.with_suffix(".png"))


if __name__ == "__main__":
    main()
