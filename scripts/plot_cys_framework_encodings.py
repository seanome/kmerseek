"""Draw what each cysteine-aware alphabet keeps of a shared cysteine framework.

Two real disulfide-rich peptides, both 36 residues, both with the six cysteines in the
same positions, from animals in different phyla:

  * O-superfamily conotoxin MiXIVA, Conus miles (cone snail, Mollusca), UniProt P06751
  * U1-theraphotoxin-Pm1a, Pterinochilus murinus (baboon spider, Arthropoda), P0DQZ2

Their shared cysteine pattern is written C-C-CC-C-C: six cysteines in four groups, with
stretches of other residues between them. Each is read three ways:

  * cys2, which writes cysteine as `c` and all 19 other residues as `x`
  * hp_pbotc_1st_ed_hpc3 (`cys_hp3`), hydrophobic `h` / polar `p` / cysteine `c`
  * hp_pbotc_1st_ed2, hydrophobic `h` / polar `p`, with cysteine counted hydrophobic

The figure is the evidence for the claim in the README's HP family section: cys_hp3
carries the cysteine positions and the hydrophobic/polar pattern of the stretches
between them at the same time, which cys2 drops and the two-class reading blurs.

Peptide sequences and the framework string come from
2024-kmerseek-analysis-cys-frameworks/tables/cys_frameworks/queries.csv, which holds
5124 mature peptides from ConoServer and Swiss-Prot.

Run from the repository root:

    python scripts/plot_cys_framework_encodings.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

import pubfig as pf  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs/images/cys_framework_encodings_conotoxin_P06751_spidertoxin_P0DQZ2"

# Mature peptides, one-letter code. Both 36 residues with cysteines at the same positions.
PEPTIDES = [
    ("Conus miles P06751", "snail", "ECREKGQGCTNTALCCPGLECEGQSQGGCAWTISGV"),
    ("P. murinus P0DQZ2", "spider", "ECRWFWGGCNNDADCCKHLECKRKWPHICLWDGTFT"),
]
FRAMEWORK = "C-C-CC-C-C"

# The partitions as src/rust/alphabets.rs defines them, residues in alphabetical order.
PBOTC_H = "ACFILMPVWY"  # hp_pbotc_1st_ed2: cysteine counted hydrophobic
PBOTC_P = "DEGHKNQRST"
HPC_H = "AFILMPVWY"  # hp_pbotc_1st_ed_hpc3: cysteine pulled out of h
HPC_P = "DEGHKNQRST"

# A realistic k for a reduced alphabet, and short enough that a 36-residue peptide has
# windows at it: kmerseek is run at k of 15 and up in hydrophobic/polar space.
KSIZE = 15

CYS_TINT = "#F3D9A4"  # a lighter tint of Okabe-Ito orange, for a cysteine column
MATCH = pf.OKABE_ITO["blue"]  # the two peptides encode to the same symbol here


def encode(sequence, classes):
    """One symbol per residue. `classes` maps a symbol to the residues it covers."""
    table = {residue: symbol for symbol, residues in classes.items() for residue in residues}
    return "".join(table[residue] for residue in sequence)


ALPHABETS = [
    (
        "cys2",
        "c = C, x = the 19 others",
        {"c": "C", "x": "ADEFGHIKLMNPQRSTVWY"},
        "Every position agrees, so cys2 cannot tell these two peptides apart at all.",
    ),
    (
        "hp_pbotc_1st_ed_hpc3 (cys_hp3)",
        f"h = {HPC_H}, p = {HPC_P}, c = C",
        {"h": HPC_H, "p": HPC_P, "c": "C"},
        "The six c symbols are the framework; the h and p between them tell the peptides apart.",
    ),
    (
        "hp_pbotc_1st_ed2",
        f"h = {PBOTC_H} (C among them), p = {PBOTC_P}",
        {"h": PBOTC_H, "p": PBOTC_P},
        "Tells them apart just as well, but the cysteines are 6 of 14 h in the snail and 6 of 17\n"
        "in the spider, so the framework cannot be read off the string.",
    ),
]


def shared_windows(first, second, ksize):
    """Windows of length `ksize` that encode identically in the two peptides."""
    windows = [
        (first[start : start + ksize], second[start : start + ksize])
        for start in range(len(first) - ksize + 1)
    ]
    return sum(a == b for a, b in windows), len(windows)


def main():
    pf.use_style()
    length = len(PEPTIDES[0][2])
    assert all(len(sequence) == length for _, _, sequence in PEPTIDES)
    cysteine_columns = [
        position
        for position in range(length)
        if all(sequence[position] == "C" for _, _, sequence in PEPTIDES)
    ]
    assert len(cysteine_columns) == 6, cysteine_columns

    rows_per_block = 5  # residues, encoded, match bars, encoded, residues
    block_gap = 4.2
    height_mm = 40 + len(ALPHABETS) * (rows_per_block + block_gap) * 4.0
    fig, ax = pf.figure(pf.TWO_COLUMN_MM, height_mm)
    ax.set_axis_off()

    # Rows above the first block hold the claim line and the legend, in reading order.
    header_rows = 7.5
    bottom = header_rows + (len(ALPHABETS) - 1) * (rows_per_block + block_gap) + rows_per_block + 3
    ax.set_xlim(-11, length + 7)
    ax.set_ylim(bottom, 0)

    mono = {"family": "monospace", "fontsize": 6.5}
    label = {"fontsize": 6, "ha": "right", "va": "center"}

    # Residue numbers every 10 positions, above the first block.
    number_row = header_rows - 0.9
    for position in range(9, length, 10):
        ax.text(position, number_row, str(position + 1), ha="center", va="center", fontsize=5)
    ax.text(-0.6, number_row, "1", ha="center", va="center", fontsize=5)
    ax.text(-10.5, number_row, "residue number", fontsize=5, ha="left", va="center")

    for block, (alphabet, classes_line, classes, note) in enumerate(ALPHABETS):
        top = header_rows + block * (rows_per_block + block_gap)
        encoded = [encode(sequence, classes) for _, _, sequence in PEPTIDES]
        agree = [a == b for a, b in zip(*encoded)]
        shared, windows = shared_windows(*encoded, KSIZE)

        # Cysteine columns, drawn first so the letters sit on top of the tint.
        for position in cysteine_columns:
            ax.add_patch(
                Rectangle(
                    (position - 0.5, top - 0.45),
                    1,
                    rows_per_block - 0.1,
                    facecolor=CYS_TINT,
                    edgecolor="none",
                    zorder=0,
                )
            )

        ax.text(
            -10.5,
            top - 1.9,
            f"{alphabet}   {classes_line}",
            fontsize=6,
            fontweight="bold",
            ha="left",
            va="center",
        )
        ax.text(
            -10.5, top - 2.8, "abc"[block], fontsize=8, fontweight="bold", ha="left", va="center"
        )

        lines = [
            (PEPTIDES[0][2], "snail residues"),
            (encoded[0], "snail encoded"),
            (None, None),
            (encoded[1], "spider encoded"),
            (PEPTIDES[1][2], "spider residues"),
        ]
        for row, (text, row_label) in enumerate(lines):
            y = top + row
            if text is None:
                for position, same in enumerate(agree):
                    if same:
                        ax.text(
                            position, y, "|", color=MATCH, ha="center", va="center", **mono
                        )
                ax.text(
                    -0.8,
                    y,
                    f"{sum(agree)}/{length} same symbol",
                    color=MATCH,
                    fontsize=5.5,
                    ha="right",
                    va="center",
                )
                continue
            for position, symbol in enumerate(text):
                ax.text(position, y, symbol, ha="center", va="center", **mono)
            ax.text(-0.8, y, row_label, **label)

        ax.text(
            length + 0.8,
            top + 1.5,
            f"{shared}/{windows} windows\nidentical at k={KSIZE}",
            fontsize=5.5,
            ha="left",
            va="center",
        )
        ax.text(-10.5, top + rows_per_block - 0.4, note, fontsize=5.5, ha="left", va="top")

    legend = [
        Patch(facecolor=CYS_TINT, edgecolor="none", label="position holding a cysteine in both peptides"),
        plt.Line2D([], [], color=MATCH, lw=1, label="the two peptides encode to the same symbol here"),
    ]
    ax.legend(
        handles=legend,
        loc="upper left",
        bbox_to_anchor=(-10.5, 2.6),
        bbox_transform=ax.transData,
        ncol=1,
        frameon=False,
        fontsize=6,
        handlelength=1.4,
    )
    ax.text(
        -10.5,
        0.2,
        f"Two toxins sharing the {FRAMEWORK} cysteine framework: cys2 reads them as the same "
        "string, cys_hp3 tells them apart and still shows\nwhere the cysteines are, and "
        "hp_pbotc_1st_ed2 tells them apart but hides the cysteines among the other "
        "hydrophobic residues.\n"
        "Conotoxin MiXIVA (Conus miles, snail, UniProt P06751) and "
        "U1-theraphotoxin-Pm1a (Pterinochilus murinus, spider, P0DQZ2), 36 residues each, "
        "from\n5124 mature peptides in ConoServer and Swiss-Prot.",
        linespacing=1.6,
        fontsize=6,
        va="top",
        ha="left",
    )

    pf.save(fig, OUT)
    print(f"wrote {OUT}.png / .pdf / .svg")
    for alphabet, _, classes, _note in ALPHABETS:
        encoded = [encode(sequence, classes) for _, _, sequence in PEPTIDES]
        shared, windows = shared_windows(*encoded, KSIZE)
        same = sum(a == b for a, b in zip(*encoded))
        print(f"{alphabet}: {same}/{len(encoded[0])} positions, {shared}/{windows} windows at k={KSIZE}")
        for line in encoded:
            print("   ", line)


if __name__ == "__main__":
    main()
