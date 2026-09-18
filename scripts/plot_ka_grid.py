import sys, math
sys.path.insert(0, __import__("os").path.dirname(__file__))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_ka_survival import load, slope_standard_error, OBSERVED, REFERENCE, FITTED, WINDOW, FLOOR, MIN_BIN_COUNT

"""Grid of calibration curves: one null per pair of rows, one database per column, from the
curves under docs/images/ka_fit_curves. `--with-database` adds the database null on top.

    python scripts/plot_ka_grid.py [--with-database] -o docs/images/ka_fit_grid.png
"""
CURVES = "docs/images/ka_fit_curves"
datasets = [("SCOPe40 (15,177 domains)", CURVES + "/scope40_{}.csv"),
            ("Swiss-Prot, 15,000 of 575,748", CURVES + "/swissprot_sample15000_{}.csv"),
            ("UniRef50, 15,000 of 38.8 M", CURVES + "/uniref50_sample15000_{}.csv")]
nulls = [("shuffled keeping dipeptides", "shuffled_dipeptide"), ("shuffled", "shuffled"), ("reversed", "reversed")]
include_database = "--with-database" in sys.argv
if include_database:
    nulls = [("database (default), grey = same queries shuffled keeping dipeptides", "database")] + nulls
out = sys.argv[sys.argv.index("-o") + 1] if "-o" in sys.argv else "ka_fit_grid.png"

fig, axes = plt.subplots(2 * len(nulls), 3, figsize=(13.5, 3.2 * 2 * len(nulls)), dpi=130)
for i, (null_label, suffix) in enumerate(nulls):
    for j, (db_label, pattern) in enumerate(datasets):
        path = pattern.format(suffix)
        meta, scores, survival, density, window, fit_density, fit_survival, ln_intercept, ref_density, ref_survival, width = load(path)
        lam, k = float(meta["slope"]), float(meta["k"])
        se = slope_standard_error(scores, density, window, lam, ln_intercept)
        xmax = 30
        for r, (observed, fitted, reference, marker, face) in enumerate([(density, fit_density, ref_density, "o", ""), (survival, fit_survival, ref_survival, "s", "none")]):
            ax = axes[2 * i + r][j]
            ax.axvspan(min(window), max(window) + width, color=WINDOW, lw=0)
            if any(reference):
                keep = [t for t, o in enumerate(reference) if o > 0 and scores[t] <= xmax]
                ax.plot([scores[t] for t in keep], [reference[t] for t in keep], "^", ms=3, mfc="none", color=REFERENCE)
            keep = [t for t, o in enumerate(observed) if o > 0 and scores[t] <= xmax]
            style = dict(mfc=face) if face else {}
            ax.plot([scores[t] for t in keep], [observed[t] for t in keep], marker, ms=3, color=OBSERVED, **style)
            keep = [t for t in range(len(scores)) if scores[t] <= xmax]
            ax.plot([scores[t] for t in keep], [fitted[t] for t in keep], "-", color=FITTED, lw=1.6)
            ax.axhline(MIN_BIN_COUNT, color=FLOOR, ls=":", lw=1)
            ax.set_yscale("log"); ax.set_ylim(0.5, 3e7); ax.set_xlim(0, xmax)
            ax.spines[["top", "right"]].set_visible(False)
            if r == 0:
                ax.set_title(f"{null_label} | {db_label}\nslope = {lam:.3f} ± {se:.3f}, K = {k:.4f}", fontsize=9, loc="left")
            if j == 0:
                ax.set_ylabel("regions in the bin" if r == 0 else "regions at or above the bin", fontsize=9)
            if r == 1 and i == len(nulls) - 1:
                ax.set_xlabel("x = λ_pair · S (nats)")
handles = [plt.Line2D([], [], marker="o", color=OBSERVED, ls="", ms=4, label="regions in each 0.5-nat bin of x (observed)"),
           plt.Line2D([], [], marker="s", color=OBSERVED, mfc="none", ls="", ms=4, label="regions at or above the bin (same regions summed)"),
           plt.Line2D([], [], color=FITTED, lw=1.6, label="fitted line, slope −(λ scale); K from the intercept"),
           plt.Rectangle((0, 0), 1, 1, color=WINDOW, label="the 8 bins the line is fitted on (≥ 30 regions each)"),
           plt.Line2D([], [], color=FLOOR, ls=":", lw=1, label="30 regions: bins below this are not fitted")]
if include_database:
    handles.insert(2, plt.Line2D([], [], marker="^", color=REFERENCE, mfc="none", ls="", ms=4, label="same queries shuffled keeping dipeptides (reference for where the database fit stops)"))
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.985 if include_database else 0.98), ncol=2, frameon=False, fontsize=9)
fig.suptitle("Calibration queries under each null (rows) against each database (columns); hp_thomas_dill2 k = 12, C = 2, X = 8, 200 queries\n"
             "x = λ_pair · S: the pair's closed-form λ times the region score S = matches − 2 × mismatches; slope 1 means the closed form holds",
             fontsize=11, x=0.02, ha="left", y=1.0)
fig.tight_layout(rect=(0, 0, 1, 0.955 if include_database else 0.945))
fig.savefig(out, bbox_inches="tight")
if out.endswith(".png"):
    fig.savefig(out[:-4] + ".svg", bbox_inches="tight")
print(out)
