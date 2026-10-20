#!/usr/bin/env python3
"""Region E-values before and after counting Pr(same class) from the region, drawn.

    cargo test --no-default-features --lib write_region_pr_same -- --ignored
    python scripts/plot_region_pr_same_nb344.py

(a) Notebook 344's worst run E-value (seanome/2024-kmerseek-analysis PR 174): a dipeptide
shuffle of DSPP against the KE repeat of TBC1D10B, 60 positions. Both rows are almost all
polar (p) in hp_lehninger2, so inside the region two positions share a class almost
always, while the two whole proteins (DSPP 15% hydrophobic, TBC1D10B 53%) give a chance of
0.48.

(b) The E-value each rule reports (`region_evalue`: the Karlin-Altschul one when there is
one, else the run E-value) for the two notebook 344 decoys and for CED-9 against BCL2, a
real match. hp_lehninger2, mismatch penalty 2, give-up margin 8, K 0.0008 and r_database
0.738 as notebook 344's index fitted, m and n of the human proteome (11,395,293 residues,
20,600 proteins). Decoys at k = 19, CED-9/BCL2 at k = 15. With SEG masking neither decoy
keeps a region.

(c) The cost on ordinary sequence: the 25 BCL-2-like test proteins against 500 decoys
(each shuffled 20 times keeping its dipeptide counts), exact search at hp_lehninger2
k = 15. Query-target pairs whose best region has a run E-value at or below each cutoff. A
calibrated E-value gives 12,500 (1 - e^(-x/500)) pairs, about 25 x while x is small.
"""

import argparse
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pubfig as pf  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

DATA = "docs/images/region_pr_same_nb344"
DECOYS_FASTA = "tests/testdata/fasta/nb344_decoys_dspp_eln_dipeptide_r2.fasta"
TARGETS_FASTA = "tests/testdata/fasta/nb344_targets_tbc1d10b_fmn2.fasta"
# 0-based region of the DSPP decoy on TBC1D10B, from region_run_evalue_holds_in_a_low_complexity_region.
DSPP_START, DSPP_END, TBC_START = 1057, 1117, 715
# hp_lehninger2's classes, from src/rust/alphabets.rs (LEHNINGER_HP).
HP = str.maketrans("AFGILMPVWYCDEHKNQRST", "hhhhhhhhhhpppppppppp")
N_PAIRS, N_TARGETS = 25 * 500, 500

# One colour per rule, the same in panels b and c.
RULES = {
    "protein": (pf.OKABE_ITO["vermillion"], "Pr(same) from the two whole proteins (before)"),
    "region": (pf.OKABE_ITO["blue"], "Pr(same) from the region (this change, the default)"),
    "protein + SEG": (pf.OKABE_ITO["bluish_green"], "whole proteins, SEG-masked seeds (--seg-mask)"),
}
CASES = [
    ("DSPP decoy x TBC1D10B", "DSPP decoy against TBC1D10B\n(chance match, 60 positions, k = 19)"),
    ("ELN decoy x FMN2", "elastin decoy against FMN2\n(chance match, 407 positions, k = 19)"),
    ("CED-9 x BCL2", "CED-9 against BCL2\n(real match, 26 positions, k = 15)"),
]


def read_tsv(name):
    with open(os.path.join(DATA, name)) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path):
    records, name = {}, None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = ""
            elif line:
                records[name] += line
    return records


def plot_residues(ax):
    decoy = read_fasta(DECOYS_FASTA)["Q9NZW4__dipeptide_r2"][DSPP_START:DSPP_END]
    target = read_fasta(TARGETS_FASTA)["sp|Q4KMP7|TB10B_HUMAN"][TBC_START : TBC_START + len(decoy)]
    rows = [
        (f"DSPP decoy {DSPP_START + 1}", decoy, f"{DSPP_END}"),
        ("  hp_lehninger2", decoy.translate(HP), ""),
        ("  same class", "".join("|" if a == b else " " for a, b in zip(decoy.translate(HP), target.translate(HP))), ""),
        ("  hp_lehninger2", target.translate(HP), ""),
        (f"TBC1D10B {TBC_START + 1}", target, f"{TBC_START + len(decoy)}"),
    ]
    # One monospaced string per row, so the label, residues and end number never overlap.
    width = max(len(label) for label, _, _ in rows)
    for i, (label, seq, end) in enumerate(rows):
        ax.text(0, 1 - i * 0.2, f"{label:<{width}}  {seq}  {end}", family="monospace", ha="left", va="center")
    n_p = (decoy.translate(HP).count("p"), target.translate(HP).count("p"))
    ax.text(
        0,
        -0.12,
        f"hp_lehninger2, k = 19: p = CDEHKNQRST (polar), h = AFGILMPVWY. {n_p[0]} of 60 decoy and "
        f"{n_p[1]} of 60 target positions are p.\n"
        "Pr(same class) from the two whole proteins: 0.477, run E-value 4.2e-10 (called significant). "
        "From the region's own positions: 0.715, run E-value 7.5 (not significant).",
        ha="left",
        va="center",
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.25, 1.1)
    ax.axis("off")


def plot_evalues(ax):
    rows = {(r["case"], r["rule"]): r for r in read_tsv("examples.tsv")}
    offsets = {"protein": -0.22, "region": 0.0, "protein + SEG": 0.22}
    no_region_x = 1e9
    for i, (case, _) in enumerate(CASES):
        ax.axhline(i, color=pf.GREY, lw=0.3, zorder=0)
        for rule, (colour, _) in RULES.items():
            r = rows[(case, rule)]
            y = i + offsets[rule]
            if r["evalue"] == "NA":
                ax.plot(no_region_x, y, "x", color=colour, ms=4, mew=1)
                continue
            e = float(r["evalue"])
            ax.plot(e, y, "o", color=colour, ms=4)
            # Labels sit on the side of the dot away from the E = 1 line.
            label = f"{e:.1e}  (Pr(same) {float(r['pr_same']):.2f})"
            if e < 1:
                ax.text(e / 12, y, label, va="center", ha="right")
            else:
                ax.text(e * 12, y, label, va="center")
    ax.axvline(1, color="black", lw=0.5, ls="--")
    ax.set_yticks(range(len(CASES)), [label for _, label in CASES])
    ax.invert_yaxis()
    ax.set_xscale("log")
    ticks = [1e-30, 1e-20, 1e-10, 1, 1e6]
    ax.set_xticks(ticks + [no_region_x], ["1e-30", "1e-20", "1e-10", "1", "1e6", "no\nregion"])
    ax.set_xlim(1e-34, 1e11)
    ax.set_xlabel("region_evalue (dashed: E = 1; smaller = more significant)\nhp_lehninger2, penalty 2, human-proteome m and n")


def plot_decoys(ax):
    rows = read_tsv("decoys_pairs_at_or_below.tsv")
    xs = [10 ** (i / 20) for i in range(-10, 61)]
    ax.plot(xs, [N_PAIRS * (1 - math.exp(-x / N_TARGETS)) for x in xs], color="black", lw=1)
    for rule in ("protein", "region"):
        colour = RULES[rule][0]
        mine = [r for r in rows if r["rule"] == rule and int(r["pairs"]) > 0]
        ax.plot([float(r["evalue_cutoff"]) for r in mine], [int(r["pairs"]) for r in mine], "o-", color=colour, ms=2.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([0.3, 1, 10, 100, 1000], ["0.3", "1", "10", "100", "1000"])
    ax.set_yticks([1, 10, 100, 1000, 10000], ["1", "10", "100", "1000", "10,000"])
    ax.set_xlim(0.3, 1000)
    ax.set_xlabel("run E-value cutoff x\n25 BCL-2-like proteins against\n500 decoys, hp_lehninger2 k = 15")
    ax.set_ylabel("decoy pairs with best region at E ≤ x\n(of 12,500)")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", default="docs/images/region_pr_same_nb344")
    args = parser.parse_args()

    pf.use_style()
    fig = pf.plt.figure(figsize=(pf.TWO_COLUMN_MM / 25.4, 120 / 25.4), layout="constrained")
    top, bottom = fig.subfigures(2, 1, height_ratios=[1, 1.9])
    a = top.subplots()
    b, c = bottom.subplots(1, 2, width_ratios=[1.6, 1])
    plot_residues(a)
    plot_evalues(b)
    plot_decoys(c)
    handles = [Line2D([], [], color=colour, marker="o", ls="", ms=4, label=label) for colour, label in RULES.values()]
    handles += [
        Line2D([], [], color="black", marker="x", ls="", ms=4, mew=1, label="no region left (SEG masked its seeds)"),
        Line2D([], [], color="black", lw=1, label="(c) calibrated E-value: 12,500 (1 - e^(-x/500))"),
    ]
    bottom.legend(handles=handles, loc="outside upper center", ncol=2, frameon=False)
    for ax, letter in zip((a, b, c), "abc"):
        pf.panel_label(ax, letter)
    pf.save(fig, args.output, formats=("png",))


if __name__ == "__main__":
    main()
