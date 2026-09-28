"""Draw the walk that grows one matched region past its exact seed, on both sides.

`extend_regions` in src/rust/search.rs walks outward from each end of a seed along the
encoded sequences: +1 for a position where the two classes agree, minus the mismatch
penalty where they differ. It stops once the running score has fallen the give-up margin
(BLAST calls this the X-drop) below its best so far, and keeps each side up to that best.
The score starts at 0 on the seed edge, so a side whose score never rises above 0 keeps
nothing, however far the walk went before giving up.

Replays that walk in Python on the JSON `kmerseek pair` writes, for the longest region,
and draws the residues, their classes, and the running score on both sides:

    kmerseek pair --query tests/testdata/fasta/bcl2.fasta --target tests/testdata/fasta/ced9.fasta \\
        --ksize 12 --alphabet hp --output scripts/testdata/bcl2_vs_ced9.hp.k12.pair.json
    python scripts/plot_xdrop_walk.py scripts/testdata/bcl2_vs_ced9.hp.k12.pair.json \\
        --output docs/images/xdrop_walk_bcl2_ced9_bh1
"""

import argparse
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

SEED_COLOR = "#9e9e9e"  # exact seed
KEPT_COLOR = "#2a9d8f"  # residues the walk adds to the region
GIVE_UP_COLOR = "#e76f51"  # the give-up line
SCORE_COLOR = "#222222"  # the running score


def walk(query_enc, target_enc, pairs, penalty, xdrop):
    """One side of the walk.

    `pairs` yields (query index, target index) stepping away from the seed. Returns
    (kept, steps): `kept` is how many positions the best-scoring extension covers, and
    `steps` holds one (query index, target index, running score) per position visited,
    including the ones after the best that the give-up margin threw away.
    """
    score = best = 0.0
    kept = 0
    steps = []
    for n, (qi, ti) in enumerate(pairs, start=1):
        score += 1.0 if query_enc[qi] == target_enc[ti] else -penalty
        steps.append((qi, ti, score))
        if score > best:
            best, kept = score, n
        elif best - score > xdrop:
            break
    return kept, steps


def walk_both_ways(pair, region, penalty=2.0, xdrop=8.0):
    """The left and right walks from one region of a `kmerseek pair` report."""
    q, t = pair["query"]["encoded"], pair["target"]["encoded"]
    qs, qe = region["query_start"], region["query_end"]
    ts, te = region["target_start"], region["target_end"]
    right = walk(q, t, zip(range(qe, len(q)), range(te, len(t))), penalty, xdrop)
    left = walk(q, t, zip(range(qs - 1, -1, -1), range(ts - 1, -1, -1)), penalty, xdrop)
    return left, right


def short_name(fasta_name):
    """`sp|P10415|BCL2_HUMAN Apoptosis regulator ...` -> `BCL2_HUMAN`."""
    token = fasta_name.split()[0]
    parts = token.split("|")
    return parts[2] if len(parts) >= 3 else token


def title_line(pair, left_kept, right_kept, penalty, xdrop):
    return (
        f"The walk grows the seed {right_kept} residues to the right and {left_kept} to the "
        f"left ({short_name(pair['query']['name'])} vs {short_name(pair['target']['name'])}, "
        f"{pair['moltype']}, k={pair['ksize']}; mismatch penalty {penalty:g}, "
        f"give-up margin {xdrop:g})"
    )


def draw_letters(ax, pair, region, left_steps, right_steps):
    """Five rows of text: query residue, query class, agree/differ, target class, target residue."""
    q_raw, t_raw = pair["query"]["sequence"], pair["target"]["sequence"]
    q_enc, t_enc = pair["query"]["encoded"], pair["target"]["encoded"]
    offset = region["target_start"] - region["query_start"]
    columns = [qi for qi, _, _ in left_steps] + list(
        range(region["query_start"], region["query_end"])
    ) + [qi for qi, _, _ in right_steps]
    for qi in columns:
        ti = qi + offset
        agree = q_enc[qi] == t_enc[ti]
        rows = [q_raw[qi], q_enc[qi], "|" if agree else "x", t_enc[ti], t_raw[ti]]
        for y, ch in zip(range(4, -1, -1), rows):
            ax.text(qi, y, ch, ha="center", va="center", family="monospace", fontsize=10)
    labels = [
        f"{short_name(pair['query']['name'])} residue",
        "its class",
        "same class?",
        "its class",
        f"{short_name(pair['target']['name'])} residue",
    ]
    ax.set_yticks(range(4, -1, -1), labels, fontsize=9)
    ax.set_ylim(-0.6, 4.6)
    ax.tick_params(axis="both", length=0)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)


def draw_side(ax, edge, kept, steps, xdrop, pair, region):
    """One walk: its running score, the give-up line under it, and the residues it kept."""
    q_enc, t_enc = pair["query"]["encoded"], pair["target"]["encoded"]
    xs = [edge] + [qi for qi, _, _ in steps]
    scores = [0.0] + [s for _, _, s in steps]
    ax.plot(xs, scores, color=SCORE_COLOR, lw=1.4, zorder=3)
    for qi, ti, s in steps:
        marker = "o" if q_enc[qi] == t_enc[ti] else "x"
        ax.plot(qi, s, marker, color=SCORE_COLOR, ms=6, zorder=4)
    best_so_far = []
    best = 0.0
    for s in scores:
        best = max(best, s)
        best_so_far.append(best)
    ax.plot(xs, [b - xdrop for b in best_so_far], "--", color=GIVE_UP_COLOR, lw=1.2, zorder=2)
    if kept:
        direction = 1 if steps[0][0] > region["query_start"] else -1
        lo, hi = sorted((edge, edge + direction * kept))
        ax.axvspan(lo, hi, color=KEPT_COLOR, alpha=0.3, lw=0, zorder=1)


def annotate(ax, left, right, region, xdrop):
    (left_kept, left_steps), (right_kept, right_steps) = left, right
    qs, qe = region["query_start"], region["query_end"]
    left_end = left_steps[-1]
    right_end = right_steps[-1]
    ax.text(left_end[0] + 0.2, 2.6, f"never above 0: keeps {left_kept}", fontsize=9)
    ax.annotate(
        f"stops {abs(left_end[2]):g} below its best",
        xy=(left_end[0], left_end[2]),
        xytext=(left_end[0] + 0.6, left_end[2] - 1.6),
        fontsize=9,
        arrowprops={"arrowstyle": "-", "color": SCORE_COLOR, "lw": 0.8},
    )
    peak_score = max(s for _, _, s in right_steps)
    peak_x = right_steps[right_kept - 1][0]
    ax.annotate(
        f"best {peak_score:+g} after {right_kept} residues: keeps {right_kept}",
        xy=(peak_x, peak_score),
        xytext=(peak_x - 1, 2.6),
        fontsize=9,
        arrowprops={"arrowstyle": "-", "color": SCORE_COLOR, "lw": 0.8},
    )
    ax.annotate(
        f"stops {peak_score - right_end[2]:g} below its best",
        xy=(right_end[0], right_end[2]),
        xytext=(right_end[0] - 6, right_end[2] - 1.6),
        fontsize=9,
        arrowprops={"arrowstyle": "-", "color": SCORE_COLOR, "lw": 0.8},
    )


def legend_handles(penalty, xdrop):
    return [
        Line2D([], [], color=SCORE_COLOR, lw=1.4, label="running score"),
        Line2D([], [], color=SCORE_COLOR, marker="o", ls="", label="same class: +1"),
        Line2D([], [], color=SCORE_COLOR, marker="x", ls="", label=f"different class: -{penalty:g}"),
        Patch(color=SEED_COLOR, alpha=0.35, label="exact seed"),
        Patch(color=KEPT_COLOR, alpha=0.3, label="residues the walk adds"),
        Line2D([], [], color=GIVE_UP_COLOR, ls="--", label=f"give-up line: best so far minus {xdrop:g}"),
        Line2D([], [], color=SCORE_COLOR, lw=0.5, alpha=0.5, label="score 0, where each walk starts"),
    ]


def plot_walk(pair, region, penalty=2.0, xdrop=8.0):
    left, right = walk_both_ways(pair, region, penalty, xdrop)
    (left_kept, left_steps), (right_kept, right_steps) = left, right
    qs, qe = region["query_start"], region["query_end"]
    fig, (letters, ax) = plt.subplots(
        2, 1, figsize=(13, 6.2), sharex=True, gridspec_kw={"height_ratios": [1.2, 2.6]}
    )
    for panel in (letters, ax):
        panel.axvspan(qs - 0.5, qe - 0.5, color=SEED_COLOR, alpha=0.35, lw=0, zorder=1)
    draw_letters(letters, pair, region, left_steps, right_steps)
    draw_side(ax, qs - 0.5, left_kept, left_steps, xdrop, pair, region)
    draw_side(ax, qe - 0.5, right_kept, right_steps, xdrop, pair, region)
    annotate(ax, left, right, region, xdrop)
    ax.text(
        (qs + qe) / 2 - 0.5,
        -11.6,
        f"seed: {short_name(pair['query']['name'])} {qs}..{qe} = "
        f"{short_name(pair['target']['name'])} {region['target_start']}..{region['target_end']}, "
        f"{region['length']} residues, exact",
        ha="center",
        fontsize=9,
    )
    ax.axhline(0, color=SCORE_COLOR, lw=0.5, alpha=0.5, zorder=1)
    ax.set_ylim(-12.5, 3.8)
    ax.set_xlim(left_steps[-1][0] - 0.8, right_steps[-1][0] + 0.8)
    ax.set_ylabel("running score")
    ax.set_xlabel(f"position in {short_name(pair['query']['name'])} (0-based, as in the CSV)")
    fig.suptitle(title_line(pair, left_kept, right_kept, penalty, xdrop), fontsize=11, y=0.995)
    fig.legend(
        handles=legend_handles(penalty, xdrop),
        loc="upper center",
        ncol=4,
        bbox_to_anchor=(0.5, 0.965),
        fontsize=9,
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("pair_json", help="JSON written by `kmerseek pair`")
    parser.add_argument("--output", required=True, help="path without extension; writes .png and .svg")
    parser.add_argument("--mismatch-penalty", type=float, default=2.0)
    parser.add_argument("--xdrop", type=float, default=8.0, help="the give-up margin")
    args = parser.parse_args()
    with open(args.pair_json) as handle:
        pair = json.load(handle)
    fig = plot_walk(pair, pair["regions"][0], args.mismatch_penalty, args.xdrop)
    for ext in ("png", "svg"):
        fig.savefig(f"{args.output}.{ext}", dpi=150 if ext == "png" else None)


if __name__ == "__main__":
    main()
