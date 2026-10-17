---
name: writing-style
description: >
  Olga's plain, technical writing standard. Use whenever writing or editing prose of any
  kind: a README, doc, PR description, notebook markdown cell, figure legend, report,
  commit message, memory note, chat reply, or code/doc comment (///, //, docstrings). Covers cutting bolded thesis-sentence openers, meta-commentary,
  em-dash chains, filler intensifiers, and narrative framing, plus the core rule: don't
  sound smart, be smart — explain complex ideas simply instead of reaching for big words.
---

# Plain, Technical Writing

Write clearly, succinctly, and technically. Do not write to sound smart.

**This is an editing standard, not just a description — apply it.** When writing new
prose or editing existing prose in a notebook markdown cell, code/doc comment, or PR
description/body, actually rewrite the text to match these rules before considering
the writing done. Don't just avoid the habits below in new text and leave old
violations sitting nearby unfixed — if you're touching a notebook, a PR, or a
comment block for another reason and it has thesis-bolding, filler words, em-dash
chains, or unexplained jargon, fix those in the same pass. Silently applying the
rule beats calling it out inline.

## Don't sound smart, be smart

Don't reach for a big word when a plain one says the same thing. Do explain genuinely
complex or technical ideas simply — a hard concept stated clearly is smarter than the
same concept left dense. The goal isn't fewer big words for their own sake, it's
whether a reader without the background can follow the point.

Words that have slipped through before and should not again: "vacuous" (say "matches
everything, so it filters nothing"), "manufactures a reading" (say "produces"),
"interpolated" (say "read off the line between X and Y"). Added 2026-09-17.
"Orthogonal" (in the plain or the mathematical sense) and "landscape" (as in probability
landscape) are fine: Olga uses both herself. Changed 2026-09-24 from her handoff notes,
which reversed the 2026-09-17 ban on "orthogonal".

If a jargon term is necessary (e.g. "anticonservative", "multiplicity", "FWER"),
explain it inline the first time it appears, or replace it with plain wording
entirely. Don't assume the reader already knows the term.

## Habits to avoid

- **Bolded thesis sentence opening every paragraph.** Bold is for labels and genuine
  emphasis, not a rhetorical topic sentence in each block.
- **Meta-commentary about the writing or the analysis itself**: "stated rather than
  left implicit", "worth stating rather than glossing over", "the honest framing is",
  "this is the headline number", "not silently". Just state the fact; the reader
  doesn't need to be told it is being stated.
- **Em-dash chains.** Three or four clauses strung together with `—`. Use periods.
- **Filler intensifiers**: actually, genuinely, really, exactly, precisely, simply,
  cleanly, outright, decisively, remotely, comfortably.
- **Narrative framing**: "lessons learned the hard way", "tells a different story",
  "the escape hatch", "bugs before they were features", "the plot thickens".
- **Defensive hedging** that argues with an imagined reviewer instead of reporting the
  result.

## No words coined in this session

If a phrase did not exist before the conversation started ("free route", "route A",
"neighbours", "structurally unreachable", "delimit"), the reader has not seen it and cannot
decode it. Say the thing itself, in the words a biologist already uses: "a different Pfam
domain on the same protein", "a call covering the whole protein". A variable name is not a
word either: `n_solo_targets` is "a target protein whose only Pfam domain is this family".
Added 2026-09-16 after one explainer needed five rounds of clarification.

## No word with two meanings

A word that already names one thing in the project cannot be borrowed for another, in
code or in a sentence. In the kmerseek repos "sweep" is the alphabet x ksize parameter
sweep. A Makefile target `sweep-dark-set` that cancelled leftover SLURM jobs, and a
handoff line saying "run the sweep first", read as launching another parameter sweep.
Grep the project for a word before naming a target, flag, function or file with it, and
before using it in a handoff. Say the action itself: cancel, delete, clean up, rebuild.
Added 2026-09-17 after the `sweep-dark-set` target had to be renamed `cancel-orphan-tasks`.

## Keep the facts intact

Keep every number, filename, gene symbol, citation, and technical qualifier intact
when editing — the problem is the prose around the facts, not the facts.

## Where this applies

READMEs, docs, PR descriptions, notebook markdown cells, reports, and code/doc
comments (`///`, `//`, docstrings) — anywhere prose explains something to a reader.
One idea per sentence. Lead with the finding, not with a characterization of the
finding.
## Keep Olga's plain names; put the field's term in parentheses once

When a review says "the field's word is X-drop" or "κ is Cohen's κ", do not replace the
plain name Olga chose ("the give-up margin", "the copy rate") with the jargon. Keep the
plain name everywhere and add the field's term once, where the thing is first defined:
"X, the give-up margin (BLAST calls this the X-drop)". A plain name that says what the
thing does beats a name the reader has to look up, and a name Olga picked is hers to keep.
Added 2026-09-20 after a review pass swapped both of those out and she asked for them back.

## Check a report or document before sending it

`scripts/check_prose.py <file.html|.md|.txt> [--max 30]` lists every sentence over the
word limit (grouped by MultiQC section id for a rendered report, tables and column-chooser
modals skipped) and counts the cut-list words. Run it on the rendered HTML, not the
builder: the report-fixes list of 2026-09-20 counted 187 long sentences that the builder's
line breaks hid. Added 2026-09-20 when the fixes file referred to "the check script in the
writing-style skill" and there was none.

## Label a fold by what is in it

A `<details>` summary or any "more" toggle names its contents in plain words: "What the
five metrics measure, and how a rank is counted". Never "Show the working", "Details",
"More" or "Under the hood". The reader decides whether to open it from the label alone.
Added 2026-09-23 after Olga struck "Show the working" from the alphabet explainer.

## Math notation

Added 2026-09-24 from Olga's handoff notes, after she called "h·h + p·p" terrible.

- Set math in real LaTeX (`$...$` / `$$...$$`) where it renders: notebook markdown cells,
  HTML pages and artifacts. In a chat window that does not render `$...$`
  (Olga, 2026-09-25), put each formula on its own line in a plain ```text``` code
  block, with the worked numbers underneath.
- Name every quantity. One symbol, one value: write `h_query`, `p_target`, not a bare `h`
  used for two proteins.
- In hydrophobic-polar alphabet math, `h` and `p` mean only the hydrophobic and polar
  class shares. Every probability is `Pr(·)`: `Pr(same)`, never `p_same`.

## Citations

Added 2026-09-24 from Olga's handoff notes; she has been burned by invented identifiers.

- Never write a DOI, PMID or URL from memory. Resolve it against the source (PubMed,
  Crossref, the publisher page) or give author, title, journal and year and say the DOI
  still needs checking.
- Chicago style, journal in italics: Karlin, S., and S. F. Altschul. "Title." *Journal
  Full Name* 87, no. 6 (1990): 2264–68. https://doi.org/...
- An annotation on a reference goes in a separate bold note under the entry, as Nature
  Reviews does, not mixed into the citation.

## Words for a biology reader

Do not call a grid entry a "cell" (in a heatmap, table or matrix): the readers are
biologists. Say "square", "entry" or name what it holds ("the chicken, k19 value").
Added 2026-09-24 from Olga's handoff notes.

## Measured values keep their column name

A name that is a kmerseek output column (`region_mean_idf`, `region_evalue`,
`region_ka_u`) means only the value the code computed against the real database. An
estimate of the same thing gets a different name, a stated source, and a check against the
measured values, and never sits in the same table unmarked. Added 2026-09-24 from Olga's
handoff notes.

## Show the residues

Whenever a sequence pair, region or window is discussed, print the real aligned residues:
both sequences, coordinates, a match line and the count ("25 of 32 identical"). For a
reduced-alphabet match, put the encoded strings underneath with their own match line.
kmerseek alignments are ungapped, so never draw or describe one with gaps. Added
2026-09-24 from Olga's handoff notes.

## No `#<number>` unless it means a PR or issue

On GitHub, `#2` in a comment, review, PR body, issue body or commit message links to PR or
issue 2 of that repo. Escaping with `\#2` or `&#35;2` does not stop it. So a rank or label
never takes a `#`: write "Reviewer 2", "step 3", "No. 2", "the top hit", not "Reviewer #2"
or "the #1 hit". When you mean a PR, write "PR #56" or `seanome/kmerseek#60`. Added
2026-10-01 after the analysis-review heading "Reviewer #2" linked to PR 2 in PR #56's body,
and a scan found "TTN is the #1 hit" in PR #83. Check every comment, PR body and commit message
for a stray `#` before posting it.
