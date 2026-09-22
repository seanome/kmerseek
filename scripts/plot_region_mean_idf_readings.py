"""Draw how region_mean_idf adds up for a region through ambiguous residues.

Topi pancreatic ribonuclease (P00659) against bovine RNase A (P61823), matched
region residues 38-102, protein20 at k=10, against the 125 UniProtKB ribonucleases
in tests/testdata (121 distinct sequences, the count the index keeps). A window
holding n residues that are B (Asp or Asn) or Z (Glu or Gln) is sketched under 2^n
readings; one of them is the k-mer bovine holds. The figure stacks the IDF of every
reading per window and shows three ways to average them.

IDF of a k-mer is ln(N / number of target proteins holding it), N = 121. The sums
reproduce the search CSV's region_tfidf for this region: 1512.11 for topi and
188.47 for goat RNase (P67926), which is topi with every B and Z resolved.

Run from the repository root:

    python scripts/plot_region_mean_idf_readings.py
"""

import gzip
import itertools
import math
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent.parent
FASTA = ROOT / "tests/testdata/fasta/ribonuclease_125_entries_uniprotkb_2026_09_17.fasta.gz"
OUT = ROOT / "docs/images/region_mean_idf_readings.png"

KSIZE = 10
# Region residues 38-102 (1-based, inclusive): window starts 38..93, 56 windows.
REGION_FIRST_WINDOW, REGION_LAST_WINDOW = 38, 93
ALTERNATIVES = {"B": "DN", "Z": "EQ", "J": "IL"}


def read_fasta(path):
    sequences = {}
    name = None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = ""
            else:
                sequences[name] += line
    return sequences


def readings(kmer):
    """Every k-mer an ambiguous k-mer stands for, in protein20."""
    return ["".join(t) for t in itertools.product(*[ALTERNATIVES.get(c, c) for c in kmer])]


def kmers(sequence):
    out = set()
    for i in range(len(sequence) - KSIZE + 1):
        out.update(readings(sequence[i : i + KSIZE]))
    return out


def main():
    sequences = read_fasta(FASTA)
    # The index keeps one copy of identical sequences, so frequencies count distinct ones.
    distinct = list(dict.fromkeys(sequences.values()))
    n_targets = len(distinct)
    frequency = Counter()
    for sequence in distinct:
        frequency.update(kmers(sequence))
    topi = next(s for n, s in sequences.items() if "RNAS1_DAMKO" in n)
    bovine = kmers(next(s for n, s in sequences.items() if "RNAS1_BOVIN" in n))

    def idf(kmer):
        return math.log(n_targets / frequency[kmer])

    windows = []
    for start in range(REGION_FIRST_WINDOW, REGION_LAST_WINDOW + 1):
        kmer = topi[start - 1 : start - 1 + KSIZE]
        windows.append([(idf(r), r in bovine) for r in readings(kmer)])
    n_windows = len(windows)
    n_readings = sum(len(w) for w in windows)
    idf_all = sum(v for w in windows for v, _ in w)
    idf_matched = sum(v for w in windows for v, matched in w if matched)

    blue, grey = "#2f6fb3", "#555555"
    fig, (ax, table) = plt.subplots(
        2, 1, figsize=(10, 6.4), gridspec_kw={"height_ratios": [3, 1.1]}, constrained_layout=True
    )
    for start, window in zip(range(REGION_FIRST_WINDOW, REGION_LAST_WINDOW + 1), windows):
        bottom = 0.0
        # Matched reading at the bottom of the stack, then the others.
        for value, matched in sorted(window, key=lambda q: not q[1]):
            if matched:
                ax.bar(start, value, bottom=bottom, width=0.8, color=blue, edgecolor=blue, linewidth=0.5)
            else:
                ax.bar(start, value, bottom=bottom, width=0.8, color="none", edgecolor=grey, linewidth=0.6)
            bottom += value
    ax.set_xlabel(
        f"window start in topi ribonuclease (aa); {n_windows} windows, "
        f"{REGION_FIRST_WINDOW}–{REGION_LAST_WINDOW}, all matched by bovine"
    )
    ax.set_ylabel("sum of IDF over the window's readings (nats)")
    ax.set_xlim(REGION_FIRST_WINDOW - 1, REGION_LAST_WINDOW + 1)
    ax.set_ylim(0, 82)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(
        handles=[
            Patch(facecolor=blue, edgecolor=blue, label="reading bovine RNase A holds (one per window)"),
            Patch(facecolor="none", edgecolor=grey, label="other reading of the same window, not in bovine"),
        ],
        loc="upper right",
        frameon=False,
        title=f"box height = IDF, ln({n_targets} / proteins holding the k-mer);\n"
        f"{math.log(n_targets):.1f} when topi alone holds it",
        title_fontsize=9,
        fontsize=9,
    )
    ax.set_title(
        "region_mean_idf adds the IDF of every reading at a position, then divides by the number of positions",
        fontsize=11,
        loc="left",
    )

    table.axis("off")
    rows = [
        ("Today: every reading, divided by windows", f"all boxes ÷ {n_windows}",
         f"{idf_all:.0f} ÷ {n_windows} = {idf_all / n_windows:.1f}"),
        ("Matched reading only, divided by windows",
         f"filled boxes ÷ {n_windows}; what goat RNase scores",
         f"{idf_matched:.0f} ÷ {n_windows} = {idf_matched / n_windows:.2f}"),
        ("Every reading, divided by readings", "all boxes ÷ number of boxes",
         f"{idf_all:.0f} ÷ {n_readings} = {idf_all / n_readings:.2f}"),
    ]
    for k, (title, detail, value) in enumerate(rows):
        y = 0.85 - k * 0.33
        table.text(0.01, y, title, fontsize=10.5, fontweight="medium", va="center", transform=table.transAxes)
        table.text(0.01, y - 0.13, detail, fontsize=9, color=grey, va="center", transform=table.transAxes)
        table.text(0.99, y - 0.06, value, fontsize=10.5, fontweight="medium", va="center", ha="right",
                   transform=table.transAxes)
        table.plot([0, 1], [y - 0.2, y - 0.2], color="#cccccc", lw=0.5, transform=table.transAxes)
    fig.savefig(OUT, dpi=170)
    print(f"{OUT}: {n_windows} windows, {n_readings} readings, IDF sum {idf_all:.2f}, matched {idf_matched:.2f}")


if __name__ == "__main__":
    main()
