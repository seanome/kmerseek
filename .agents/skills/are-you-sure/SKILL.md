---
name: are-you-sure
description: Adversarial "are you sure?" pass before calling any technical work done. Use after writing or editing any script, pipeline process, Nextflow config, notebook code cell, plotting code, Rust, or SQL/polars query; before handing Olga any shell command; before reporting what a command, squeue/sacct, a log or a file on Sherlock shows; and before any git commit, push, merge, branch change or PR. Also use when Olga asks "are you sure?", "did you check?", or "does that actually work?". Every coding skill points here.
---

# Are you sure? The adversarial pass on code, commands, Sherlock and git

Added 2026-09-30 at Olga's request. On 2026-09-29 one "Are you sure?" on the experiments page
found 7 claims that a 570-number audit had passed. Code has the same problem: it runs, it
prints something plausible, and the bug is in what it did not do. This pass is the step where
you try to prove your own code wrong before Olga has to ask.

Run it every time, after the tests pass and before you say the work is done. Passing tests is
the start of this pass, not the end.

## 1. Say what the code claims

Write down, for yourself, one line per claim the code makes. "This reads every combo",
"this retries on OOM", "this filters to the window", "these ranks are 1..n with no ties".
Each claim is something to attack in step 2.

## 2. Attack each claim with a real input

For each claim, find the input that would break it, then run the code on that input or
read the code path that handles it. Reasoning "it should be fine" does not count. The
questions that have found real bugs in these repos:

- **Empty, one, many.** Zero rows, one file, one match. A glob that matches one file gives a
  bare path, not a list. An empty result can be written to a cache and look like a real zero.
- **Missing values.** What happens to nulls after a left join? polars `min_horizontal` and
  `max_horizontal` skip nulls, which turned every false positive into IoU 1.0.
- **Ties and order.** Does anything depend on row order? Sorting on score alone leaves ties to
  chance. A lazy `with_row_index` is not stable. Print `n_tied` for any rank.
- **Silent success.** Does a failure look like success? `cmd | tail` hides the exit status.
  `|| true` stored 178 OOM kills as empty results. BioMart serves an outage as HTTP 200. A run
  that searched 1_000 of 19_696 queries looks like real negatives.
- **The case the guard was written for.** Trigger the error path, not only the happy path. A
  Nextflow task killed by the scheduler has exit status `Integer.MAX_VALUE`, so a retry on
  `128..143` never fired on the OOM it was written for.
- **Scale.** Does it still fit in memory and time on the real input, not the test subset?
  `scan_csv` on a `.csv.zst` inflates the whole file in RAM. For a cluster launch: has
  every task's first memory AND time ask been checked against the earlier runs' failures,
  including the run being resumed? Set each task's ask from the
  measured history, not from a retry. A task that would need a retry is a bug in the launch.
- **Stale input.** Is it reading the file you think? A default path pointing at an old
  snapshot, a cached result from before your fix, an image tag that was never pushed.
- **Names you did not resolve.** Every function, flag, tag, URL and column name checked
  against the source or the pinned docs, not recalled.
- **Merges.** After merging branches that grew the same file, check for duplicate `def`s.

## 3. Check the output against something you did not compute

Compare at least one number the code produced against an independent count: `wc -l`, a
second query written a different way, a row picked by hand, the number in the paper. If the
code and the check agree, say which check. If you could not find one, say so.

## 4. Read the diff as a reviewer who wants to reject it

`git diff` the whole change and read it top to bottom as if someone else wrote it. Look for
debugging leftovers, a hard-coded path or ksize, a changed default that other callers rely
on, a comment that no longer matches the code, and work outside the task.

## 5. Report what the pass found

In the message that says the work is done, include one short line per claim: what you tried
and what happened. Name anything you could not test ("did not run on Sherlock", "no
scheduler kill to test against"). Never write "verified" for something only reasoned about.

## 6. Commands handed to Olga

Added 2026-09-30 when Olga extended this pass to all command-line work. Before a command goes
in a `bash` block for her:

- Run it yourself first when it is safe to (read-only, or on a scratch copy). If you could not,
  say so next to the command.
- Say where it runs: her Mac, a Sherlock login node, or inside a Sherlock job. Start it with the
  full `cd /absolute/path &&`. On Sherlock give the bare command, no `ssh sherlock` wrapper.
- Check it against her shell: `cp`, `mv`, `rm` are aliased to `-i` and hang or no-op, so use
  `/bin/cp -f` etc.; `ls` is a recursive `lsd`, use `/bin/ls`; in zsh `$VAR:r`, `:h`, `:t` are
  modifiers, write `${VAR}`; start any `gh` command with `GH_PAGER=cat`.
- Sherlock's git is 1.8.3.1: no `git -C`, `switch` or `restore`.
- Any code the command pulls or runs is pushed first. Check with `git log origin/<branch> -1`
  that the commit she will get is the one you tested.
- Prefer the `make` target over the raw command when one exists. Restate the whole command every
  time; never "run the same as before".

## 7. Reading command output, local or on Sherlock

Before you report what a command showed ("the array finished", "no errors", "3 files"):

- Look at the raw output once, unfiltered, before summarising it.
- Check that the field you grep for is in the output at all. `squeue -o "%i %t %M %R" | grep 244`
  can never match a job name, because `%j` is not in the format. An empty result from a broken
  filter looks the same as an empty queue.
- Count units, not lines. A pending array range `45021171_[68-99]` is one `sacct` row and 32
  tasks; use `squeue -r` to expand it.
- Check the output is complete: a pager, `head`, a truncated log or a `--limit` default can cut
  it. Output that ends mid-list is cut off.
- Check you read the right thing: the right clone (Sherlock has a second clone), the right run
  directory, the newest trace (Nextflow truncates its trace at run start), today's file.
- When Olga's screen says something different from your reading, your reading is the suspect.

## 8. Git: commits, pushes, merges, branches

Before any commit:

- `git status` and `git diff --staged`, read in full. Stage named paths only, never `git add -A`
  or `git add .`: Olga edits the same worktrees.
- Run the tests or the check without a pipe, or with `set -o pipefail`: `pytest | tail && git commit`
  commits even when the tests fail.
- Confirm the branch (`git branch --show-current`) and the worktree you are in.
- The commit message says what changed and why, matches the diff, and ends with the attribution
  line.

Before a push, merge or branch change:

- `git log --oneline origin/<branch>..HEAD` to see exactly what will go up.
- Never read another branch's files with `git checkout <branch> -- <path>`: it writes the index.
  Use `git show <branch>:<path>` or a worktree.
- After a merge, check for duplicate top-level `def`s in files both sides changed.
- After pushing, read back what landed (`git log origin/<branch> -1`, or the PR body via
  `gh api ... --jq .body`), not what you meant to send.
- Never switch branches under a running Nextflow head: `bin/` is read at task time.

If the pass finds a bug, fix it and run the pass again on the fix. If the same kind of bug
turns up twice, add it to the list in step 2 with the date and what it cost.
