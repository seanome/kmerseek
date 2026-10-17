---
name: clear-figures
description: Rules and a pre-publish checklist for any figure, diagram, plot, SVG, matplotlib panel, or interactive explainer. Use whenever drawing anything visual for Olga (an inline diagram, an HTML page, a notebook figure, a figure for a paper or talk), and whenever explaining a concept, because every explanation gets a picture. Also use when the user says a figure is unclear, asks what a mark or word means, or says "I don't understand this section".
---

# Clear figures, no coined words

Olga is a visual learner. She reads the picture before the caption and the legend
before the text. A figure that needs its caption to be understood has failed.

Added 2026-09-16 after one explainer app went through five versions because grey meant
"rest of the protein" on one row and "another domain" on the next, and the prose used
words invented during the session ("free route", "neighbours", "delimit") that the reader
had no way to decode. Both mistakes had been made earlier in the same session.

Fonts, sizes, line weights, colours and export (Nature Biotechnology specs) are in the
`publication-figures` skill. This skill is about what the marks mean; run both.

## Encoding rules

1. **One mark per meaning, one meaning per mark.** Two different things never share a
   colour, shape, or line style. If you find yourself explaining in a caption that "grey
   here means X but grey there means Y", the figure is wrong. Change a mark.
2. **Every mark is in a legend, and the legend comes before the marks** (above or to the
   left, in reading order). If a mark is not in the legend, remove it or add it.
3. **Use the field's own drawing conventions.** Proteins are a thin line for the backbone
   with boxes for domains (the Pfam / InterPro convention). Never a filled bar for the
   protein and another filled bar for a domain. Alignments and calls are bars under the
   protein line. Trees are trees. Genomes are lines with arrows for genes.
4. **Real units on every axis and every mark**: aa, bp, Mya, %, n. A bar without a number
   is decoration.
5. **Colour carries meaning, never sequence.** Do not cycle through a palette because
   there are several items. Group by category, and use at most three hues plus a neutral.
   Semantic colours (pass / fail, present / absent) are separate from category colours.
6. **Pass/fail, present/absent, and other binary states get a shape or position cue as
   well as a colour**, so the figure survives greyscale and colour-blindness.
7. **Labels do not collide and do not leave the canvas.** Compute the text width
   (~7 px per character at 12 px) against the space available before placing it. Long
   labels move to a legend or below the row; they do not shrink.
7a. **Never write "(log scale)" in an axis label.** It reads as if the values were
    transformed. Label the ticks with the real values and what they mean (scaled
    "1 (all)", "2 (1/2)", "5 (1/5)", "10 (1/10)"), and if the spacing is by ratio, say so
    in plain words once above the figure: "the step from 1 to 2 takes as much room as the
    step from 5 to 10". Use plain tick numbers (20, 50, 100), not 10^2. Added 2026-10-01
    after Olga read "scaled (log scale)" in notebook 270 and could not tell what it meant.

Rules 7c to 7h come from Olga's PR review comments, collected 2026-10-01 (PR numbers are
seanome/2024-kmerseek-analysis unless marked kmerseek).
In 2024-kmerseek-analysis, 7c-7f are built into `notebooks/figure_utils.py`
(seanome/2024-kmerseek-analysis#91, draft as of 2026-10-01): `alphabet_axis(ax, alphabets)` groups alphabets by size
and puts grid lines at the ticks, `MISSING_STYLE` is the missing-value dot, and
`finish_figure` (also re-exported by `mhc_region_utils` and `hp_conservation_utils`)
measures the top of the plot and puts the TOOLS line right above it. Use them instead of
redoing these by hand.

7c. **Grid lines start at the tick marks.** When a panel has one row or column per
    category (an alphabet, a species), draw its guide line from the axis tick straight
    across the panel, so every mark can be read back to its label. A grid offset from the
    ticks, or a grid with no ticks, makes the reader guess the row. PR #46 asked for this
    on two figures of notebook 241.
7d. **Group a category axis by what the categories are, with a gap between groups.**
    Alphabets go in size groups: 2-3 letters, 4-8, 12-18, then protein20 alone, with blank
    space between groups. Not alphabetical, not in the order the data came in. Same for
    species (by clade) and tools (kmerseek arms, then the other tools). PR #46.
7e. **"No value" gets a mark that looks like nothing else.** A missing result (no
    E-value fit, no hit, not run) is a small dot `.` or a grey ×, never a circle when
    the data are circles, and never a hollow version of a data marker. PR #46: the
    "no E-value fit" circle read as another data point.
7f. **No empty band between the header text and the plot.** The tool line, hypothesis
    and conclusion at the top of a figure sit directly above the legend and axes. Olga
    flagged the gap on one figure of notebook 241 as "this figure and honestly all of
    them" (PR #46). `tight_layout` leaves a fixed top margin and says nothing (see the
    project memory `figures_tight_layout_silently_refuses`): measure `ax.get_position()`
    and the text's bounding box, and close the gap by hand.
7g. **Show the parameters the figure depends on.** A figure of a hydrophobic-polar (HP)
    alphabet lists which residues are hydrophobic and which polar, since that differs by
    alphabet (kmerseek PR #58). When several alphabets share a core, show the shared part
    once in code font (`AFILMV` hydrophobic, `DEHKNQRST` polar) and list only the
    residues that move (kmerseek PR #43).
7h. **Real proteins at real sizes.** Examples use proteins Olga knows (BCL2/Ced9,
    human-mouse GENCODE genes) at the k-mer sizes actually run: HP k of 19 or more,
    never k=9 ("completely unreasonable", kmerseek PR #62; kmerseek PR #38).

Rules 7i to 7r come from Olga's chat messages, August to October 2026 (collected
2026-10-01), including the MultiQC figure spec she pasted on 2026-09-03.

7i. **Numbers go on an x-y plot, not in bars.** "Anytime the categories are actually
    numeric, show me an x vs y plot instead of bars" (2026-09-03): k, divergence time,
    disorder, length, class count. Plot each item at its measured value, not at the
    midpoint of a bin: "Can this figure not be by disorder bins but by real disorder of the
    region?" (2026-09-04).
7j. **Never connect categories with a line.** A line says the values in between exist.
    Alphabet, species, feature type and tool get points (the 2026-09-03 spec).
7k. **Show the whole population, mark the pick.** Never plot only the winner ("which
    alphabet won"). Draw every alphabet-ksize pair as a dot or heatmap square and highlight
    the chosen one, so a lead inside the noise is visible as one (2026-09-03 spec).
7l. **Show the spread, and name the summary.** A mean over nine species gets the nine
    values beside it (points, range or box), and the label says mean, median or total:
    "What is 'pooled'? Mean, median? Can we see a range/boxplot ... for all species?"
    (2026-09-03).
7m. **Few points are points.** No violin or density curve over n below ~30: a curve over 19
    values invents a shape. Draw the points and label the ones the text names (2026-09-01).
7n. **Colour scale fitted to the data.** A colorbar runs over the data's own range, not a
    fixed 0-1 that turns the panel into one flat colour. Sequential colormap for a
    quantity with no natural middle (Fmax, recall); a diverging one only around a real
    zero, such as a difference (2026-09-01). The light end of a sequential map must still
    show against a white background: start it at a mid tint (`cmap(np.linspace(0.25, 1))`)
    or give each point a thin dark edge. Added 2026-10-01 after the checklist, tried on
    notebook 241's P66/CD47 rank figure, found ranks near 10,000 drawn almost white.
7o. **Use colour.** Olga asked for it outright ("remember to use colors", 2026-09-10;
    "What happened to the alternating colors?", 2026-09-29). Greyscale-only panels and
    dropped colour bands are regressions. Rules 5 and 6 still apply.
7p. **Name the tool, and for kmerseek the alphabet and k, on every panel** (2026-09-10),
    and say which protein set is the query and which the target database (2026-09-17).
7q. **Show examples a biologist knows.** Next to a summary of how far off the calls are,
    show a few named genes (2026-09-10: "Could we see some examples of 'famous' genes?").
    SHOW the result, don't describe it (2026-09-10, about synteny).
7r. **Say how to read the figure and what to decide from it.** "What is this figure really
    telling me? How do I interpret it? What decisions do I make from it?" (2026-09-25). The
    conclusion line answers the last question.

## Words in and around figures

7b. **A mark that another mark covers is a missing mark.** Two series that nearly
    coincide (a metric and the control it turns out to equal) will overplot, and the one
    drawn last wins. Draw the reference or control *first*, as a wide pale band, and the
    series on top of it thinner and dashed, so both stay readable. Check by naming every
    series in the legend and finding each one in the rendered image before publishing.
    Added 2026-09-22, after region_tfidf was invisible under the region_length control
    it correlated with at rho 0.967 -- the whole point of the panel, drawn and then hidden.

8. **No words coined in this session.** If a phrase did not exist before you started
   ("free route", "route A", "neighbours", "structurally unreachable", "delimit"), the
   reader has not seen it either. Say the thing itself: "a different Pfam domain on the
   same protein", "a call covering the whole protein".
9. **Name things by what a biologist calls them**, not by what the code calls them.
   `n_solo_targets` is "a target protein whose only Pfam domain is this family".
10. **Define a necessary technical term the first time it appears, in the same sentence.**
    IoU: "the overlap between the call and the true domain divided by their union".
11. **The one fact the figure rests on goes in the title or first line**, not in the
    third paragraph. If the whole picture depends on "a call's interval is copied from
    the target protein it aligned to", say that first.
12. **Say which direction is good.** A hypothesis or conclusion line that reports a
    number says whether higher is better or worse, and what value means "no signal". PR
    #46: "Is 'expected chance' bigger than 0 good or bad?" A term like "whether a lambda
    can exist at all" gets the plain meaning in the same line, or goes.
13. **Symbols match the explainer, and each one is named.** Use the letters the published
    explainer uses (`u` for the match probability, not `a`; kmerseek PR #68), and never
    leave a bare `r`, `t` or `hi` on a figure or in its code (kmerseek PR #54, #68). Plain
    words for terms like "censored" and "decoy", defined where they first appear
    (kmerseek PR #54, #71).
14. **A figure stands on its own, and lives somewhere.** Someone without the PR context
    must be able to read it; link the explainer page for the method rather than assume it
    (kmerseek PR #73). Its file name says what it shows and from which data
    (`xdrop_walk_bcl2_ced9_bh1.png`, kmerseek PR #67; the source, Swiss-Prot, UniProt or
    SCOPe40, in the name, PR #52). A doc, README or notebook links it: never a bare file
    in `docs/images/` (kmerseek PR #67).
15. **A claim about a sequence or a number gets a figure, not only text.** "It's hard to
    trust just pure text" (kmerseek PR #69); "Why doesn't this grow to the left? I need
    to see a figure" (kmerseek PR #67). A section with only a table gets a figure too (PR
    #46, sections 1 and 8 of notebook 241).

## Pre-publish checklist, run every time

Before saving a figure or publishing any page or diagram that holds one:

- [ ] List every colour / shape used and the one meaning of each. Any duplicate meaning?
      Any mark missing from the legend?
- [ ] Is the legend placed before the marks in reading order?
- [ ] Are proteins drawn as line + boxes? Are units on every axis and label?
- [ ] Read every label and every sentence of caption aloud. Any word that did not exist
      before this session? Any variable name where a plain phrase belongs?
- [ ] Does every text label fit its space (width check), including the longest dynamic
      value the label can take?
- [ ] Look at the rendered output once (Read the PNG, or the widget preview). If a
      preview is impossible, say so to the user in one line.
- [ ] Find every legend entry in the rendered image. A series you cannot point to is
      covered by another one: redraw, do not explain it in the caption.
- [ ] Does any annotation or label sit on top of a curve, bar or point? Move it to
      empty space with a leader line.
- [ ] Would the figure make sense to a biologist who has not read the notebook? If the
      answer needs a caption, fix the figure, not the caption.
- [ ] Do the grid lines start at the tick marks? Is a category axis grouped (alphabet
      size, clade, tool family) with gaps between groups?
- [ ] Does a missing value use a mark no data point uses?
- [ ] Is there a white band between the header text and the legend or axes? Measure it.
- [ ] Are the alphabet's residue classes, k and data source on the figure or in its title?
      Is k realistic for that alphabet?
- [ ] Do the hypothesis and conclusion lines say which direction is good?
- [ ] Does the file name say what and from which data, and does a doc or notebook link it?
- [ ] Numeric x drawn as an x-y plot at measured values? No line joining categories?
- [ ] Whole population shown with the pick marked, spread shown next to any mean, and the
      summary named (mean, median, total, n)?
- [ ] Colorbar fitted to the data range, sequential unless there is a real zero? Colour used?
- [ ] Tool named on the panel (kmerseek: alphabet and k), query and target named?

## When the user says a figure is unclear

Do not explain the existing figure. Ask which mark or word broke, if it is not already
stated, then change the encoding so that mark can only mean one thing, and republish once.
Every clarification request in a session is a rule to add above, not a one-off patch.

## Sticky headers in interactive pages

A chart inside a horizontally scrollable box (`overflow-x:auto`) is painted on its own
layer in WebKit (Safari, iOS, and desktop apps built on WebKit) and can show through a sticky header even
when the header has the higher z-index. Give the content sections `position:relative;
z-index:0`, give each chart box `isolation:isolate`, and make the header sticky only
above ~820 px: on a phone a sticky control bar covers most of the screen. Added
2026-09-23 after a chart bar showed as a purple box on top of the alphabet selector.

## Interactive pages

From Olga's feedback on the reference-results pages, 2026-09-30:
- No wide empty margin; a side panel of alphabets gets a draggable divider so it can be
  narrowed, and the whole list fits without scrolling to the bottom.
- Hovering an alphabet name shows its residue classes.
- Picking an alphabet and k is a small heatmap or plot that shows which ones find the
  region, clicked to open the result, not two dropdowns.
- Every notebook, PR and commit mentioned is a link.
- A heatmap says what its numbers and its sort order are, in the legend above it.

## The plotting code: are you sure?

Run the `are-you-sure` pass on the code that draws the figure (added 2026-09-30),
alongside the checklist above. Also:
- Count the plotted points or bars and match the count to the table the cell prints.
- Pick one mark by hand and trace it back to its row in the data.
- Check nothing was dropped silently: nulls, values outside the axis limits, a category the
  colour map did not list.
