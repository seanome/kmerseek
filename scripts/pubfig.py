"""Helpers for Nature-family figures. Pair with nature.mplstyle (same folder).

Copy both files into the repo (2024-kmerseek-analysis: notebooks/) so a figure can be
rebuilt from committed code; do not import them from ~/.claude.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation

MM = 1 / 25.4  # inches per millimetre
ONE_COLUMN_MM = 89
TWO_COLUMN_MM = 183
MAX_HEIGHT_MM = 170  # leaves room for the legend on the page

STYLE = Path(__file__).with_name("nature.mplstyle")

# Okabe-Ito, the colour-blind-safe set Nature's figure guide lists.
OKABE_ITO = {
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "bluish_green": "#009E73",
    "reddish_purple": "#CC79A7",
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "yellow": "#F0E442",
    "black": "#000000",
}
GREY = "#999999"  # for "everything else" / context marks


def use_style():
    plt.style.use(STYLE)


def figure(width_mm=ONE_COLUMN_MM, height_mm=60, **subplots_kw):
    """Figure at its printed size, with constrained layout."""
    if height_mm > MAX_HEIGHT_MM:
        raise ValueError(f"height {height_mm} mm is over Nature's {MAX_HEIGHT_MM} mm")
    subplots_kw.setdefault("layout", "constrained")
    return plt.subplots(figsize=(width_mm * MM, height_mm * MM), **subplots_kw)


def panel_label(ax, letter, dx_pt=-18, dy_pt=4):
    """Bold lowercase 8 pt panel letter, a fixed distance up-left of the axes corner."""
    offset = ScaledTranslation(dx_pt / 72, dy_pt / 72, ax.figure.dpi_scale_trans)
    ax.text(
        0,
        1,
        letter,
        transform=ax.transAxes + offset,
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def shared_legend(fig, ncol=None, **kw):
    """One legend above all panels, de-duplicated across axes."""
    handles, labels = {}, []
    for ax in fig.axes:
        for h, lab in zip(*ax.get_legend_handles_labels()):
            if lab not in handles:
                handles[lab] = h
                labels.append(lab)
    return fig.legend(
        [handles[lab] for lab in labels],
        labels,
        loc="outside upper center",
        ncol=ncol or len(labels),
        **kw,
    )


def save(fig, path_stem, formats=("pdf", "svg", "png")):
    """Vector PDF/SVG for the journal, 450 dpi PNG for notebooks and slides."""
    path_stem = Path(path_stem)
    path_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        fig.savefig(path_stem.with_suffix(f".{ext}"))
