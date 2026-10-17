#!/usr/bin/env python3
"""List the sentences over a word limit and the cut-list words in a report or document.

    check_prose.py report.html            # a MultiQC report: grouped by section id
    check_prose.py README.md notes.txt    # any text or markdown
    check_prose.py --max 25 report.html   # a stricter limit (default 30 words)

Prints one line per hit: file, section (for HTML), word count, the sentence. Ends with
the cut-list word counts. Exit status is 0 either way; this is a reading aid, not a gate.
"""
import argparse
import html
import re
import sys
from pathlib import Path

CUT_WORDS = ["deliberately", "honest", "honestly", "manufacture", "manufactures", "simply",
             "actually", "exactly", "precisely", "truly", "delineate", "delineation",
             "the former", "the latter", "genuinely", "legible", "vacuous",
             "interpolated", "really", "cleanly", "outright", "decisively", "remotely",
             "comfortably", "worth noting", "it is worth", "note that", "importantly"]

# MultiQC chrome that is not prose.
CHROME = re.compile(r"^(AI Summary|Provider:.*|Chat with Seqera AI|Summarize .*|Copy .*|"
                    r"More details…|Export.*|Created with MultiQC|Configure columns|"
                    r"Sort by highlight|Scatter plot|Violin plot|Show (All|None)|Sort|"
                    r"Visible|Group|Column|Description|ID|Scale|\|\||Close|Table|"
                    r"Expand table|Showing .*|.*Uncheck the tick box.*)$")


def sections_of_html(text: str) -> list[tuple[str, str]]:
    text = re.sub(r"<script.*?</script>", "", text, flags=re.S)
    text = re.sub(r"<style.*?</style>", "", text, flags=re.S)
    # Tables are numbers, not prose, and MultiQC's column-chooser modals repeat every
    # column description; both would swamp the list.
    text = re.sub(r"<table.*?</table>", " ", text, flags=re.S)
    text = re.sub(r'<div class="modal.*?</div>\s*</div>\s*</div>', " ", text, flags=re.S)
    text = re.sub(r'<div class="mqc-section [^"]*" id="mqc-section-wrapper-([\w-]+)">',
                  r"\n\n@@SECTION \1\n", text)
    # Block ends become sentence ends so two <li>s are not one sentence.
    text = re.sub(r"</(p|li|h[1-6]|td|th|div)>", ". ", text)
    plain = html.unescape(re.sub(r"<[^>]+>", " ", text))
    parts = re.split(r"@@SECTION (\S+)\n", plain)
    out = [("head", parts[0])]
    for i in range(1, len(parts), 2):
        out.append((parts[i], parts[i + 1]))
    return out


def sentences(block: str):
    block = "\n".join(l for l in block.splitlines() if not CHROME.match(l.strip()))
    block = re.sub(r"\s+", " ", block)
    for s in re.split(r"(?<=[.!?])\s+(?=[A-Z\"'(])", block):
        s = s.strip(" .")
        words = s.split()
        # A run of bare numbers and dots is a table that escaped the tag filter.
        if words and sum(w in (".", "..") for w in words) < 3:
            yield s


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("paths", nargs="+", type=Path)
    ap.add_argument("--max", type=int, default=30, help="word limit (default 30)")
    args = ap.parse_args()
    counts = {w: 0 for w in CUT_WORDS}
    n_long = 0
    for path in args.paths:
        text = path.read_text(encoding="utf-8", errors="replace")
        blocks = (sections_of_html(text) if path.suffix.lower() in (".html", ".htm")
                  else [("", text)])
        for sid, block in blocks:
            for s in sentences(block):
                n = len(s.split())
                if n > args.max:
                    n_long += 1
                    where = f"{path.name}:{sid}" if sid else path.name
                    print(f"{where}\t{n} words\t{s}")
                low = " " + s.lower() + " "
                for w in CUT_WORDS:
                    counts[w] += len(re.findall(r"\b" + re.escape(w) + r"\b", low))
    print(f"\n{n_long} sentences over {args.max} words", file=sys.stderr)
    hits = {w: c for w, c in counts.items() if c}
    print("cut-list words: " + (", ".join(f"{w} {c}" for w, c in sorted(hits.items(), key=lambda t: -t[1])) or "none"),
          file=sys.stderr)


if __name__ == "__main__":
    main()
