---
name: publication-figures
description: Typography, size, line, colour and export settings for every matplotlib figure Olga will put in a paper, talk, PR or notebook, following Nature Biotechnology's figure specifications (Arial/Helvetica, 5-7 pt, 89/183 mm). Use before writing any plotting code, before any savefig, when asked for a "publication-quality", "paper", "print" or "Nature" figure, and whenever Olga comments on fonts, font size, line weight, or a figure looking unpolished. What the marks mean is the clear-figures skill; this skill is how they look.
---

# Publication figures (Nature Biotechnology style)

Added 2026-10-02. Olga had asked about fonts and sizes many times. Two open PRs show why:
PR 74's print figure came out in DejaVu Sans, and PR 101's had Arial turned into outlines,
which Nature does not accept because the text cannot be edited. One style file now sets
all of it.

Meaning rules (one mark per meaning, legend first, proteins as line + boxes) are in
`clear-figures`. Run both.

## Use the style file, do not set fonts by hand

The skill folder has two files:

- `nature.mplstyle`: every rcParam below.
- `pubfig.py`: `figure()`, `panel_label()`, `shared_legend()`, `save()`, `OKABE_ITO`, `GREY`.

In a repo, copy both into the shared code folder (kmerseek: `scripts/`; 2024-kmerseek-analysis: `notebooks/`)
and commit them, so the figure can be rebuilt from the repo. Do not import from
this skill folder. If a copy already exists in the repo, use it and do not make a second one.

```python
import pubfig as pf
pf.use_style()
fig, axs = pf.figure(pf.ONE_COLUMN_MM, 60, ncols=2)   # size in mm, as printed
...
pf.shared_legend(fig)                                  # one legend above all panels
for ax, letter in zip(axs, "ab"):
    pf.panel_label(ax, letter)
pf.save(fig, "figures/<what>_<data source>")           # .pdf .svg .png
```

Do not pass `fontsize=` or `linewidth=` in plotting calls except for the panel letter and
a deliberate emphasis line. Every hand-set size is a place the style can drift.

## The numbers (Nature figure guide, read 2026-10-02)

Source: research-figure-guide.nature.com, "Preparing figures: our specifications" and
"Building and exporting figure panels".

| What | Value | Where it comes from |
|---|---|---|
| Font | Arial, else Helvetica. Linux/Sherlock: Liberation Sans (same letter widths as Arial) | Nature |
| Text size | 5-7 pt at printed size. Axis labels 7, ticks and legend 6 | Nature range; split is ours |
| Panel letters | 8 pt bold, lowercase, upright: **a**, **b**, **c**; top left of each panel | Nature |
| Sequences | Courier (monospaced), one-letter code, 50-100 residues per line | Nature |
| Width | 89 mm (one column) or 183 mm (two columns) | Nature |
| Height | at most 170 mm | Nature |
| Lines | axes and ticks 0.5 pt, data lines 1 pt, grid 0.3 pt light grey | Nature gives no number; common choice |
| Colour | RGB; Okabe-Ito set (below); no coloured text, use a key instead | Nature |
| Text contrast | black or white text on a high-contrast background (4.5:1) | Nature |
| File | vector PDF or EPS, text kept as text, fonts embedded as TrueType (type 42) | Nature |
| Images | at least 300 dpi, 450 dpi or more preferred | Nature |
| Axes | axis lines and tick marks always; units in parentheses: "Length (aa)" | Nature |
| Avoid | patterns, drop shadows, 3-D, decoration, titles inside the figure | Nature |

Make the figure at its printed size. A figure drawn at 10 inches and shrunk to 89 mm has
2 pt text and hairlines. `pf.figure()` takes mm and refuses heights over 170 mm.

## Colour

Okabe-Ito, in the order the style cycles it: blue `#0072B2`, vermillion `#D55E00`,
bluish green `#009E73`, reddish purple `#CC79A7`, orange `#E69F00`, sky blue `#56B4E9`,
yellow `#F0E442`, black. Grey `#999999` is for context (the rest of the population,
the background). Yellow and sky blue are light: use them for fills, not thin lines.

Give a colour a meaning and keep it across all figures of the paper: if kmerseek is blue
in Figure 1, it is blue in Figure 4. Write the mapping down once in the repo
(a dict in `pubfig.py` or the notebook utils) and look colours up from it.

Continuous values: `viridis` or `cividis`; a diverging map (`RdBu_r`) only around a real
zero. Never `jet` or `rainbow`.

## Look

- Open axes: left and bottom spines only. No box, no background colour.
- No title on a paper panel. The panel's point goes in the figure legend (the caption).
  A notebook figure keeps its header line (tool, hypothesis, conclusion, from
  `finish_figure`) but the exported paper version drops it.
- Gridlines: Nature says avoid them. Keep the rule from `clear-figures` 7c, a line from
  each category tick when the reader must trace a row back to its label, and draw it at
  0.3 pt light grey. No other grid.
- Legend: no frame, above the panels (`shared_legend`) or in empty space inside one.
  A legend wider than its panel pushes the next panel away; use the shared legend.
- Markers 3 pt; scatter points with many overlaps get `s=4-8`, `alpha`, or no edge.
- Error bars say what they are (s.d., s.e.m., 95% CI) and n, in the axis label or legend.
- Plain tick numbers (20, 50, 100), never 10^2 (see `clear-figures` 7a).
- Space panels tightly and in alphabetical order, left to right, top to bottom.

## Before calling a figure done

- [ ] `pf.use_style()` called before the figure was made; no stray `fontsize=` in calls.
- [ ] Width is 89 or 183 mm; height at most 170 mm.
- [ ] No text below 5 pt. Check any label you set by hand, and annotation text.
- [ ] PDF has Arial embedded as TrueType, not outlines, not DejaVu:
      `grep -a -o "/BaseFont /[A-Za-z+-]*\|/FontFile[23]" fig.pdf | sort -u`
      must show `ArialMT` and `/FontFile2`. `DejaVu` means Arial was not found.
- [ ] No warnings about a missing font when saving (`findfont: Font family ... not found`).
- [ ] Look at the PNG: nothing overlaps, no panel pushed aside by its legend, panel letters
      line up.
- [ ] Colours come from the paper-wide mapping; the same thing has the same colour in
      every figure.
- [ ] Then the `clear-figures` checklist.

## Sherlock and Linux

Arial is not installed there; the style falls back to Liberation Sans, which has the
same letter widths. If neither is present, matplotlib uses DejaVu Sans and the PDF check
above catches it. Make the final paper PDFs on the Mac, or install `fonts-liberation`
and clear the font cache (`rm ~/.cache/matplotlib/fontlist-*.json`).
