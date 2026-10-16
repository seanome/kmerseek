"""Helpers for Nature-family figures. Pair with nature.mplstyle (same folder).

Copy both files into the repo (2024-kmerseek-analysis: notebooks/) so a figure can be
rebuilt from committed code; do not import them from ~/.claude.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
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
    """Vector PDF/SVG for the journal, 450 dpi PNG for notebooks and slides.

    The PNG is stored as 8-bit with a 256-colour palette, which keeps flat-colour figures
    a fraction of the size of a full-colour PNG. A figure with a continuous colour bar
    stays in full colour (see to_256_colours).
    """
    path_stem = Path(path_stem)
    path_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        path = path_stem.with_suffix(f".{ext}")
        fig.savefig(path)
        if ext == "png":
            to_256_colours(path)


# A matplotlib colour bar has 256 shades. On the 336 PNGs in 2024-kmerseek-analysis
# (2026-10-13), the 15 with a colour bar had 251-256 shades in one run and every other
# figure at most 29.
MIN_COLOUR_BAR_SHADES = 100


def colour_bar_shades(image, max_step=12, width=8):
    """Most shades in one straight run of pixels that looks like a colour bar.

    In the run each pixel is at most `max_step` levels from the one before it, and is
    exactly the colour of the next `width` pixels across the run, as in a bar filled
    across its width. The soft edge of a curve or a marker fails the second condition,
    so line and scatter plots score near 0. Same test as the repo's scripts/shrink_png.py.
    """
    a = np.asarray(image.convert("RGB")).astype(np.int16)
    best = 0
    for arr in (a, a.transpose(1, 0, 2)):  # bars running down, then across
        step = np.abs(np.diff(arr, axis=0)).max(axis=2)
        in_bar = step <= max_step
        for k in range(1, width + 1):
            same = np.zeros_like(in_bar)
            same[:, :-k] = (arr[1:, :-k] == arr[1:, k:]).all(axis=2)
            in_bar &= same
        shades = np.zeros(in_bar.shape[1], dtype=np.int32)
        for row_in_bar, row_step in zip(in_bar, step):
            shades = np.where(row_in_bar, shades + (row_step > 0), 0)
            best = max(best, int(shades.max()))
    return best


def to_256_colours(path):
    """Rewrite a PNG in place as 8-bit with a palette of at most 256 colours.

    A figure with a continuous colour bar is left in full colour: at 256 colours the bar
    shows visible steps (Olga, 2026-10-13).
    """
    from PIL import Image

    with Image.open(path) as image:
        if colour_bar_shades(image) >= MIN_COLOUR_BAR_SHADES:
            return
        dpi = image.info.get("dpi")
        palette = image.convert("RGB").quantize(colors=256, dither=Image.Dither.NONE)
    # quantize() drops the 450 dpi tag; without it the PNG opens at 72 dpi, 6x too large.
    palette.save(path, optimize=True, **({"dpi": dpi} if dpi else {}))
