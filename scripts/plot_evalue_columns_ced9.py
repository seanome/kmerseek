#!/usr/bin/env python3
"""Which E-value a CSV row reports: the run E-value against the Karlin-Altschul E-value
for every region CED9 shares with the 300 human BCL2-search proteins, hp_lehninger2 k=19.

    F=tests/testdata/fasta/uniprotkb_BCL2_AND_model_organism_9606_2025_02_06.fasta.gz
    kmerseek index --input $F --output idx.db --ksize 19 --alphabet hp_lehninger2
    kmerseek search --query tests/testdata/fasta/ced9.fasta --target idx.db \
        --output ced9_k19_extended.csv --ksize 19 --alphabet hp_lehninger2 \
        --min-shared-kmers 0 --max-pvalue 1 --extend-mismatch-penalty 2
    python scripts/plot_evalue_columns_ced9.py ced9_k19_extended.csv \
        -o docs/images/evalue_columns_ced9_bcl2_human_k19

Without --extend-mismatch-penalty every row's region_evalue is its region_run_evalue
(the x axis); with it, every row whose pair has a positive lambda reports
region_ka_evalue (the y axis).
"""

import argparse
import csv

import matplotlib

matplotlib.use("Agg")
import pubfig as pf

LOW, HIGH = 0.5, 200
LABELLED = {"BCL2_HUMAN", "CUL1_HUMAN", "AKT1_HUMAN"}


def read_rows(path):
    """Rows that have both E-values. region_ka_evalue is empty when the pair has no
    positive lambda, and region_run_evalue when the region could not be measured."""
    rows, skipped = [], 0
    for row in csv.DictReader(open(path)):
        if not row["region_run_evalue"] or not row["region_ka_evalue"]:
            skipped += 1
            continue
        rows.append(
            {
                "target": row["target_name"].split("|")[2].split()[0],
                "run": float(row["region_run_evalue"]),
                "ka": float(row["region_ka_evalue"]),
                "mismatches": int(row["region_n_mismatches"]),
            }
        )
    print(f"{len(rows)} rows drawn, {skipped} skipped for an empty E-value")
    return rows


def scatter(ax, rows):
    exact = [r for r in rows if r["mismatches"] == 0]
    grown = [r for r in rows if r["mismatches"] > 0]
    ax.scatter(
        [r["run"] for r in exact],
        [r["ka"] for r in exact],
        s=10,
        marker="o",
        color=pf.OKABE_ITO["blue"],
        label=f"no mismatch in region (n={len(exact)})",
    )
    ax.scatter(
        [r["run"] for r in grown],
        [r["ka"] for r in grown],
        s=12,
        marker="^",
        color=pf.OKABE_ITO["vermillion"],
        label=f"region grew past mismatches (n={len(grown)})",
    )


def guides(ax, low, high):
    ax.plot([low, high], [low, high], color=pf.GREY, linewidth=0.5, zorder=0)
    ax.axhline(1, color=pf.GREY, linewidth=0.5, linestyle=":", zorder=0)
    ax.axvline(1, color=pf.GREY, linewidth=0.5, linestyle=":", zorder=0)
    ax.annotate(
        "same value",
        (high / 2.5, high / 2.5),
        xytext=(3, -3),
        textcoords="offset points",
        ha="left",
        va="top",
        fontsize=6,
        rotation=45,
        rotation_mode="anchor",
    )
    ax.text(low, 1, " E = 1", ha="left", va="bottom", fontsize=6)


def label_targets(ax, rows):
    for r in sorted(rows, key=lambda r: r["ka"]):
        if r["target"] in LABELLED:
            ax.annotate(
                r["target"],
                (r["run"], r["ka"]),
                xytext=(5, 2),
                textcoords="offset points",
                fontsize=6,
            )
            LABELLED.discard(r["target"])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("csv")
    parser.add_argument("-o", "--output", required=True, help="path without extension")
    args = parser.parse_args()
    rows = read_rows(args.csv)
    pf.use_style()
    fig, ax = pf.figure(pf.ONE_COLUMN_MM, 80)
    ax.set_xscale("log")
    ax.set_yscale("log")
    guides(ax, LOW, HIGH)
    scatter(ax, rows)
    label_targets(ax, rows)
    ax.set_xlim(LOW, HIGH)
    ax.set_ylim(LOW, HIGH)
    ax.set_xlabel("region_run_evalue\n(longest exact run only; region_evalue without extension)")
    ax.set_ylabel("region_ka_evalue\n(whole extended region; region_evalue with extension)")
    pf.shared_legend(
        fig, ncol=1, title="CED9 vs 300 human BCL2-search proteins, hp_lehninger2 k=19"
    )
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    pf.save(fig, args.output)


if __name__ == "__main__":
    main()
