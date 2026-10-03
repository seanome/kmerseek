---
name: small-prs
description: Google's engineering practices for PR size, PR and commit descriptions, code review, and answering review comments (google.github.io/eng-practices). Use before opening any PR, when asked to split a PR or branch, when writing a PR description or commit message, when reviewing a PR, when writing or answering review comments, and when a PR grows past a few hundred lines. Olga sent the four pages on 2026-09-20 and said "remember this".
---

# Small PRs, and how to review and answer reviews

Source: https://google.github.io/eng-practices/ (small-cls, cl-descriptions,
handling-comments, standard, looking-for). Google says CL; here it is a PR.

## What "small" is

One self-contained change: one part of a feature, not the whole feature. The code still
works for users and developers after it merges. Its tests are in the same PR. A reviewer
needs nothing outside the PR, its description, the existing code, or a PR they already
reviewed. Not so small that its point is lost: a new API ships with one use of it.

Numbers: 100 lines is a reasonable PR. 1000 is usually too large. Files count too: 200
lines in one file may be fine, spread across 50 files it is not. When in doubt, cut
smaller. Reviewers rarely complain about a PR being too small.

Large is fine only for: deleting whole files (counts as one line), or output of an
automatic refactoring tool that is fully trusted.

## Why

Small PRs are reviewed faster and more thoroughly, carry fewer bugs, waste less if
rejected, merge with fewer conflicts, are easier to design well, block less, and roll back
more simply. A reviewer may reject a PR for size alone.

## How to split

- Stack: open PR 1, then start PR 2 on PR 1's branch right away. Never wait for review
  before starting the next piece.
- By files: groups that need different reviewers (a schema change and the code using it).
- Horizontally: a shared type or stub between layers lets each layer be its own PR.
- Vertically: each sub-feature is its own full-stack PR, sharing the common bits.
- Refactors go in their own PR, before the feature or fix. Moving or renaming a thing is
  never in the same PR as changing its behaviour. A local variable rename can ride along.
- Tests that only add coverage to existing code, or test helpers, can go first as their
  own PR, so a later refactor is checked against them.
- Plan the split before writing the code, not after.
- Do not break the build between stacked PRs: every PR in the stack must pass on its own.
- If it truly cannot be small (rare), ask the reviewer in advance and expect a long review.

## PR descriptions and commit messages

A description says what changed and why: the context the author had, and any decision
the code does not show. It is searched years later by someone with a faint memory of it,
so put the searchable facts in the description, not only in the code.

First line: one short imperative sentence saying specifically what this PR does, then a
blank line. "Delete the FizzBuzz RPC and replace it with the new system." It must stand
alone in a history listing. Body: the problem, why this approach, its shortcomings, and
context such as issue numbers, benchmark numbers, or design notes, with enough of the
content inline that a dead link does not lose it. Bad: "Fix bug", "Fix build", "Phase 1",
"Add convenience functions", "Moving code from A to B". Re-read the description before
merging; a PR changes during review. Tags are optional; keep them out of the first line
unless short.

## Reviewing (the standard)

Approve once the PR definitely improves the code's health, even if not perfect. There is
no perfect code, only better code. Facts and data beat opinion; the style guide beats
taste; design questions are principles, not preferences; if nothing else applies, match
the existing code. Prefix polish with `Nit:`. Read every line assigned, in context, not
just the diff window. Check: design, does it do what was meant and is that good for users,
no more complex than needed, no code for a future need, tests that fail when the code
breaks, clear names, comments that say why not what, docs updated when build/test/usage
changed. Say what was done well, not only what was wrong.

## Answering review comments

Not personal. First ask: what is the reviewer trying to tell me? If they did not
understand the code, fix the code first, then a code comment, and only last a reply in the
review tool. If disagreeing, explain the trade-off and ask which trade-off they weigh
differently, never "no". Seek consensus, then escalate rather than let a PR sit.

## A PR that changes a test's expected value

Added 2026-09-28 after PR #62 changed `14 matched regions` to `13` in a test with a one-line
code comment as the only reason; Olga asked why, then asked for a picture, then asked how
to make this happen every time.

Every changed expected value in a test (a count, a score, a row number, a fixture) gets,
before the PR is called done:

1. **The cause, found by running both builds.** Build the base branch and the PR, run the
   same command on each, and list the rows that differ. Never explain a changed number
   from reasoning alone. Check that the base build reproduces the old number exactly;
   if it does not, the old value was already stale and this PR is not the cause.
2. **A figure of the actual case**: the real residues, positions and k-mers of the pair
   whose value changed, before and after, drawn with a small script committed under
   `scripts/` and the image under `docs/images/`. Run the `clear-figures` checklist.
3. **A PR comment** (banner first line) with the figure embedded from a
   `raw.githubusercontent.com/<owner>/<repo>/<sha>/...` link, the before/after rows as a
   `text` block, the aligned residues with match counts, and the arithmetic (14 − 2 + 1 = 13).
4. **One line in the PR description** naming each test whose expected value changed and
   the old and new value.
5. **Everything else built from the same output**: docs pages, README figures, and every
   sentence that quotes one of the changed numbers (`git grep` the old number). Regenerate
   only files that differ; say which files came out identical and were left alone.

## Before opening the PR: are you sure?

Run the `are-you-sure` pass on the whole diff before opening or updating a PR (added
2026-09-30). Put its result in the PR body's testing section: one line per claim, what was
run and what it showed, and what was not tested. When reviewing someone else's PR, apply
step 2 of that skill to each change: find the input that breaks it before approving.
