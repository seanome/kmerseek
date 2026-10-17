---
name: rustacean-review
description: Expert Rustacean code review. Trigger whenever the user asks to review, critique, audit, or "look over" Rust code — a file, a diff, a PR, a branch, a crate, or a single function — or asks "is this idiomatic", "is this the Rust way", "why is this slow", "am I cloning too much", "can this be faster", or wants a pre-commit or pre-merge check on Rust. Also trigger before merging or shipping Rust, when wrapping up a Rust feature, and when the user mentions clippy, borrow-checker workarounds, allocations, clone(), lifetimes, unsafe, or MSRV/Cargo.toml concerns. Enforces these habits — verify every API against the pinned docs instead of recalling it from memory, audit each clone() and allocation, collapse redundant passes over the same data, apply Rust API Guidelines naming and error conventions, require SAFETY and why-comments, and strip conversational artifacts ("as discussed", "v2", "fixed above") out of code and comments. Also enforces what Olga keeps asking for by hand on kmerseek PRs — a figure for every test about sequences, real examples at real k-mer sizes (BCL2/Ced9 drawn residue by residue), a doc comment on every public item and field with each symbol defined, cited sources instead of magic numbers, one topic per PR, and no backward-compatibility code for unreleased index formats.
---

# Rustacean Review

Review Rust the way a senior maintainer reviews a PR from a colleague they respect: read
carefully, verify claims, flag what actually matters, and say when something is fine.

The two failure modes to avoid are **flattery** (approving code you didn't check) and
**noise** (twenty nitpicks that bury the one real bug). Both waste the author's time.

---

## Phase 0 — Ground yourself before you have opinions

Skipping this phase is what produces reviews full of confident nonsense: suggestions that
don't compile on their toolchain, APIs that don't exist, "fixes" for behavior that was
deliberate.

**1. Read `Cargo.toml` and `Cargo.lock` first.**

```bash
cat Cargo.toml                 # edition, rust-version (MSRV), features, dep versions
rustc --version && cargo --version
grep -A2 'name = "<dep>"' Cargo.lock   # the version actually in use
```

Edition and MSRV gate everything. `let ... else`, `if let` chains, async fn in traits,
`Option::is_none_or` — each has a stabilization version. A suggestion that needs 1.85 is
useless in a crate pinned to 1.75. If `rust-version` is absent, ask rather than assume the
latest.

**2. Never recall an API signature — look it up.**

Standard library and popular crates change every six weeks. What you remember may be
outdated, renamed, or from a different major version. Check the version actually in use:

- Vendored source is the ground truth: `~/.cargo/registry/src/*/<crate>-<version>/src/`
- `cargo doc -p <crate> --no-deps` then read the generated HTML
- `rustup doc --std` for offline std docs matching the installed toolchain
- `https://docs.rs/<crate>/<exact-version>/<crate>/` — always pin the version in the URL
- `https://doc.rust-lang.org/std/` and the release notes for recent stabilizations

If you cannot verify an API, do not cite it. Say "there may be a `try_` variant here worth
checking" rather than inventing `Vec::try_reserve_exact_unchecked`.

**3. Run the tools before you hand-review.**

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -W clippy::pedantic -W clippy::nursery
cargo test
cargo build --release 2>&1 | grep -i warning
```

Clippy catches `redundant_clone`, `needless_collect`, `clone_on_copy`, `needless_range_loop`
and hundreds more. Do not spend your review budget re-deriving what a linter already found —
run it, report what it says, then spend your attention on what a linter *cannot* see:
whether the design is right, whether the invariants hold, whether the comments are true.

If the tools cannot run (no toolchain, partial snippet), say so explicitly in the review and
mark findings as unverified.

**4. Read enough context to know the invariants.**

Before calling a `clone()` unnecessary, find out whether the original is used afterward.
Before calling an `unwrap()` a bug, find out whether a caller already guarantees the value.
Read the callers, the tests, and the type definitions. If you still cannot tell, that is a
finding in itself — "this invariant is not stated anywhere" is legitimate review feedback.

---

## Phase 1 — Correctness first

Nothing else matters if the code is wrong. Look for:

- **Panics on the library path** — `unwrap`, `expect`, `panic!`, `[i]` indexing, slicing,
  integer division, `unreachable!`. Each is a claim about an invariant. Is it documented?
- **Silent arithmetic** — `as` casts that truncate or wrap sign, `+`/`*` that can overflow
  in release (where overflow checks are off by default). Prefer `checked_`, `saturating_`,
  `try_from`.
- **Error handling that destroys information** — `.ok()`, `let _ =`, `map_err(|_| ...)`,
  a `?` that converts a specific error into a stringly-typed one.
- **`unsafe` without a `// SAFETY:` comment** stating which invariant makes it sound.
  Suggest `cargo miri test` for anything non-trivial.
- **Trait-contract violations** — `Hash`/`Eq` disagreeing, `Ord` inconsistent with `PartialOrd`,
  floats as map keys, code that depends on `HashMap` iteration order.
- **Concurrency** — a `MutexGuard` held across `.await`, blocking I/O inside async, `Rc`
  where `Arc` is needed, lock ordering that can deadlock, `AtomicUsize` with `Relaxed`
  ordering where the author meant to synchronize other data.

`correctness-and-style.md` has the fuller catalog with examples.

---

## Phase 2 — The clone and allocation audit

The instruction is "not too many clones," but the useful version is more precise: **every
`clone()` should be a deliberate choice, and the reviewer should be able to name which kind
it is.** Sort each one into a bucket:

| Bucket | Example | Verdict |
|---|---|---|
| Cheap and correct | `Arc::clone(&shared)`, `Rc::clone` | Fine. Idiomatic. Prefer the explicit `Arc::clone(&x)` form so readers see it's a refcount bump, not a deep copy. Say it's fine and move on. |
| Load-bearing | The owned value genuinely escapes — stored, sent to a thread, returned | Fine. Leave it. |
| Borrow-checker appeasement | Cloned to end a borrow early, or because a `&mut` overlapped a `&` | Real finding. The fix is usually restructuring: split the borrow, use indices, `std::mem::take`, or reorder statements — not a deeper clone. |
| Reflexive | `clone()` on a `Copy` type, cloning a whole struct to read one field, cloning a key to do a lookup, `String` where `&str` would do | Real finding. Cheap to fix. |

Allocation patterns worth flagging, roughly in order of how often they matter:

- **Allocating inside a hot loop** where a scratch buffer could be hoisted and `.clear()`ed
  each iteration. This is the single highest-value finding in sequence-processing code.
- `Vec::new()` + push in a loop with a known length → `Vec::with_capacity(n)`.
- `format!` where `write!` into an existing `String`/`Vec<u8>` works.
- `&Vec<T>` / `&String` parameters → `&[T]` / `&str`. Returning `&Vec<T>` → return `&[T]`.
- `.to_vec()` / `.collect()` into a temporary that is immediately only iterated → drop it.
- `Cow<'_, str>` for the common case where the value is usually borrowed and rarely modified.
- `String` keys cloned for `HashMap` lookup — `get` accepts `&str` via `Borrow`.

**Do not turn this into an ideology.** A clone in `main()`, in a test, in setup code, or in
anything that runs once is not worth a comment. Allocation matters where it repeats. If you
claim something is a performance problem, say what makes you think so (loop nesting, input
size) and recommend a benchmark rather than asserting a number you did not measure.

`allocation-and-iteration.md` has before/after pairs for each pattern.

---

## Phase 3 — The iteration audit

The target here is **redundant passes over the same data** — but be precise about what that
means, because the most common reviewer error is flagging something that isn't one.

**Not a redundant pass:** `v.iter().filter(..).map(..).sum()`. Iterator adapters are lazy
and fuse into a single traversal. Chaining is idiomatic and costs nothing.

**Actually a redundant pass:**

- Separate statements that each walk the same collection — `let a = v.iter().filter(..).count();`
  followed by `let b = v.iter().filter(..).count();`. Collapse into one loop or `fold`.
- A mid-chain `.collect::<Vec<_>>()` that exists only so the chain can continue.
- `.iter().any(..)` or `.contains(..)` *inside* a loop over another collection — that's
  O(n·m). Build a `HashSet` once, or sort and binary-search.
- `map.get(k)` followed by `map.insert(k, ..)` — two hashes of the same key. Use the
  `entry` API.
- `.len()`-indexed loops (`for i in 0..v.len() { v[i] }`) — pays a bounds check per access
  and reads worse than iterating. Clippy's `needless_range_loop`.
- `.remove(i)` inside a loop — O(n²). Use `retain`, or `swap_remove` if order is free.
- Recomputing a loop-invariant expression each iteration.
- `.chars().count()` for a length that only needs `.len()` when the data is known ASCII —
  and conversely, `.len()` used as a character count when it isn't.
- `sort()` + `dedup()` versus `HashSet`: which is right depends on whether order matters
  and whether the data is already nearly sorted. Ask, don't assert.
- **Rayon**: `par_iter()` on a short collection or trivial per-item work is a pessimization.
  Check the work per item before approving or suggesting parallelism.

---

## Phase 4 — Idiom and API shape

Judge against the [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/), not
against personal taste. The highest-value items:

- **Naming conventions carry meaning**: `as_*` is a cheap borrowed view, `to_*` is expensive
  and owned, `into_*` consumes `self`. Getting these backwards misleads every caller.
- **Take generic, return concrete**: `impl AsRef<Path>`, `&str`, `&[T]` in argument position;
  concrete types out.
- **Error types**: `thiserror` (or hand-written enums) for libraries so callers can match;
  `anyhow` for binaries. A public library API returning `Box<dyn Error>` or `String` throws
  away the caller's ability to handle anything.
- **Derive the standard traits** (`Debug`, `Clone`, `PartialEq`, and `Copy`/`Default`/`Hash`
  where they fit). Missing `Debug` on a public type is a real papercut.
- **Newtypes over primitive obsession** — a bare `usize` for a k-mer size is a bug waiting
  to be swapped with another `usize` parameter.
- Visibility: default to `pub(crate)`; `pub` is a promise.
- `#[must_use]` on constructors and pure transforms; `#[non_exhaustive]` on public enums.
- `#[inline]` only with a stated reason — it is not free and the compiler is usually right.

---

## Phase 5 — Comments, docs, and the no-archaeology rule

**Comments explain why, not what.** `// increment the counter` above `i += 1` is noise that
will eventually become a lie. The comments worth writing:

- `// SAFETY:` — mandatory above every `unsafe` block, naming the upheld invariant.
- `// INVARIANT:` — a struct field relationship the type system can't express.
- Why a non-obvious approach was chosen over the obvious one, and what breaks if it changes.
- The source of a magic number: a paper, a spec section, a measured threshold.

Doc comments on public items: one-line summary, then `# Errors` if it returns `Result`,
`# Panics` if it can panic, and a doctest example when the usage isn't self-evident.

### No conversational archaeology

Code must read as though written in one sitting by someone who has never seen the chat that
produced it. The next reader is a stranger — possibly the author in six months — and a
comment referencing a conversation is meaningless to them and starts rotting immediately.
Git already holds the history.

Flag and strip every one of these, in both existing code and anything you write:

| Artifact | Example |
|---|---|
| Conversational reference | `// as discussed`, `// per your suggestion`, `// you were right about this` |
| Iteration marker | `// v2`, `// new approach`, `// updated version`, `// second attempt` |
| Change narration | `// fixed the bug from before`, `// this replaces the old function`, `// changed from HashMap` |
| Versioned identifiers | `fn parse_v2`, `struct ConfigNew`, `process_data_fixed`, `old_handler` |
| Commented-out code kept as a record | any block of `//`'d-out previous implementation |
| Changelog headers in source files | `// 2026-03-01: switched to rayon` |
| Attribution to the assistant | `// TODO(claude)`, `// generated by AI` |
| Backward compatibility for formats nobody has used yet | `// older spellings`, aliases for names this branch just renamed, a refusal to open old indices, tests named `reduced_*` that pin the PR against itself |

The replacement is almost always: delete it, or rewrite it as a timeless statement of what
the code does and why. `// v2: now uses a HashSet because the Vec scan was O(n²)` becomes
`// HashSet rather than Vec: membership is checked once per input record, so linear scans
would dominate.`

Backward compatibility is part of this rule for kmerseek. Olga on PR #53: "I don't care about
the old indices ... I care about making the future better." The index stores the kmerseek
version; code reads that field instead of guessing from old layouts. The one exception is
the sourmash names `hp`, `protein` and `dayhoff`, which keep parsing.

---

## Phase 6 — Show it: figures, real examples, documentation

These are the comments Olga writes most often on kmerseek PRs, by hand, after the code
already compiles and passes tests. `what-olga-flags.md` has each one in her words with the PR
it came from. A missing item here goes under **Should fix**, not Consider. A PR missing the
figure or the example for its main claim is not mergeable, however clean the Rust is.

**1. A figure for every test about sequences.** "I'd really love figures for all these
tests. It's hard to trust just pure text" (#69). "Why doesn't this grow to the left? ... I
need to see a figure" (#67). For each new or changed test that asserts a region, an
extension, a score, a fit or a k-mer count:

- Is there a figure of that exact case, posted on the PR or under `docs/images/`?
- Is the file named for what it shows, e.g. `docs/images/xdrop_walk_bcl2_ced9_bh1.png`?
- Does a doc page, the README, or the doc comment of the code link to it? A bare image file
  nobody links to is a finding.
- Does it make sense to someone without the PR context, and link the explainer page when
  it shows E-value quantities (#73)?
- A script that writes HTML or a plot: is an example output on the PR (#58)?

If the figure is missing, make it (run the `clear-figures` checklist first) and attach it to
the review rather than only asking for it. When test values change, follow the `small-prs`
procedure for changed expected values (run both builds, draw the case, post it on the PR).

**2. Real examples at real parameters.**

- Sequences come from real proteins and the shared fixtures (`src/rust/tests/test_fixtures.rs`
  FASTA paths, or an existing `const` such as `BCL2_FRAGMENT_1_30`), not a new inline string
  (#54). Grep for an existing constant before accepting a new one.
- A region claim is drawn residue by residue: protein letters, encoded letters, match
  bars, both proteins, with coordinates (#50, #54):
  ```
  Ced9 pr: …RTVGNAQTD**QCPMSYGRLIGLISFGGFV**AAKMMESVE…
  Ced9 hp: …pphhphppp**pphhphhphhhhhphhhhh**hhphhpphp…
                     |||||||||||||||||||
  BCL2 hp: …hhphhpphh**pphhphhphhhhhphhhhh**phpphppph…
  BCL2 pr: …FATVVEELF**RDGVNWGRIVAFFEFGGVM**CVESVNREM…
  ```
- k-mer sizes are ones she would run: protein about 5, Dayhoff about 7, HP alphabets 15 or
  more and usually 19+. "k=9 for hp space ... is completely unreasonable" (#62).
- Tests cover every alphabet and every pair, not a sample or a sliding pair (#43, #56). Build
  the list from the alphabet registry and assert its length.
- Any alphabet shown names its hydrophobic and polar residues, since they differ (#58).
- Amino acids listed in a test, constant or table are in a standard order: alphabetical by
  one-letter code (`ACDEFGHIKLMNPQRSTVWY`), or grouped by the alphabet's classes with the
  grouping named. Never an ad hoc order. "It's not even alphabetical or like grouped by
  anything ... can you pick a normal order?" (#112, added 2026-10-04).
- A new behaviour comes with the case that shows why it matters, e.g. "examples where the
  whole protein passes but the region doesn't" (#38), and says how it changes downstream
  numbers such as the region score and E-value (#57).

**3. Documentation on everything public.**

- Every `pub` item, struct field, enum variant and CLI flag has a doc comment. "Make sure
  ALL properties have a doc comment. Nothing should be empty" (#54).
- A tunable parameter's doc comment gives the default and why (#54: mismatch penalty 2,
  X-drop 8).
- Every symbol and term is defined where it first appears. She has had to ask what `r`,
  `t`, `hi`, "censored", "non-modal", "multiplicity", "decoy" and "closed form" mean. A
  one-letter name is allowed only when it is the symbol the explainer page uses, and the doc
  comment says so (#68: `match_prob`, "called u in the explainer").
- Formulas are written out in the doc comment and the PR text: "lambda is the sum of
  database abundances of the k-mers from the region, divided by the region size" (#38);
  "the slope is K and intercept is lambda" (#54).
- A number that is not a p-value is not named like one (#38).
- A non-obvious line gets a why-comment (#40, #38: "what does saturating_sub help with
  here?"). A variable name says what the value is for, not that someone "wanted" it (#40).
- New `docs/*.md` pages are reachable from seanome.github.io/kmerseek or the README (#75),
  and `CHANGELOG.md` has an entry (#50).

**4. No magic numbers in tests.** "I really need to see the original citation to trust it.
Give the originals and no magic numbers!" and "This applies to ALL the tests here" (#54).
Tests still assert exact literals, and each literal says where it came from: the paper and
table, the hand calculation written out, or the command that produced it. A number used in
more than one test becomes a named constant (#54: "Can the 242 be a variable"). Expected
values pasted in from a run with no derivation "feel like cheating" (#54).

---

## Phase 7 — Scope, structure and defaults

- **One topic per PR.** "Why are all these RocksDB changes in this PR?" (#54). Changes
  outside the PR's title go to their own PR; see the `small-prs` skill.
- **Code lives in the module it belongs to.** Logic in `main.rs` that the CLI does not call
  directly moves to its module (#54).
- **No repeated fields or near-duplicate functions.** Several structs carrying the same
  fields should share one nested struct (#54: "can't they be inherited or derived? This is
  hard to maintain"). A struct with one field probably wants its siblings in it (#68). Two
  functions that load the same thing become one (#53, #63).
- **Use a library before writing one**, and weigh how heavy the dependency is against what
  it is used for (#54, #58).
- **The index is the source of settings.** A flag the user must repeat at search time, when
  the index already knows the value, is a finding (#68, #71). The default is the option she
  would choose (#70: shuffled dipeptide decoys).
- **kmerseek domain rules.** Self-matches are always skipped (#63). Ambiguous residues are
  counted by windows, not by expanded k-mers (#74).
- **No new pass over the whole index for a summary number** that could be counted while
  the index is built (#36, #37). Shared state written from parallel code is checked for
  races (#48). "How much time does this add?" is answered with a release-build timing, not
  a guess (#37).

---

## Output format

Use this structure. It puts the author's attention where it belongs.

```markdown
## Verdict
[2–3 sentences: is this mergeable, and what's the one thing that matters most.]

## Tooling
[What you ran and what it said. Or: what you couldn't run and why.]

## Blocking
[Correctness bugs, unsound unsafe, panics on reachable paths. Often empty — good.]

## Should fix
[Real problems that will bite: wrong API shape, hot-loop allocations, misleading comments.]

## Figures, examples and docs
[Phase 6 checklist, one line per item: present (with link), missing (with what to add), or
made by this review (attach the figure). Never leave this section out.]

## Consider
[Judgment calls and taste. The author may reasonably decline.]

## Verified clean
[What you checked and found genuinely good. Name it specifically.]
```

Each finding gets:

1. **Location** — `src/sketch.rs:142`
2. **What** — one sentence, no hedging
3. **Why it matters** — the concrete consequence, not "this is bad practice"
4. **The fix** — a minimal diff or snippet that compiles on *their* toolchain
5. **Confidence** — if the finding depends on an invariant you couldn't verify, say which

Example finding:

> **`src/sketch.rs:142` — `String` allocated per k-mer inside the scan loop.**
> `encode(&kmer.to_string())` builds a fresh `String` for every window; for a 10 kb
> sequence at k=26 that's ~10,000 allocations per record. `encode` only reads the bytes,
> so it can take `&[u8]` and borrow the window directly.
> ```rust
> // before
> for kmer in seq.windows(k) { out.push(encode(&kmer.to_string())); }
> // after
> for kmer in seq.windows(k) { out.push(encode(kmer)); }
> ```
> Confidence: high on the allocation; the speedup depends on how `encode` dominates the
> profile — worth a criterion run before and after.

---

## Record the review on the PR

When the code under review is a GitHub PR, or a branch that has an open PR, post the finished
review as a comment on that PR. Olga asked for this on 2026-09-28 so she can see which PRs have
had a review, and the weekly PR status review reads the marker to fill its "rustacean-review"
column. This is a standing request: post without asking first. Post only the review itself, and
only on the PR that was reviewed.

1. Find the PR and the exact commit reviewed:
   ```bash
   gh api repos/<owner>/<repo>/pulls/<n> --jq '.head.sha'
   ```
   If the working tree has commits not yet pushed, or the head SHA on GitHub differs from what
   you read, say so in the comment. The marker must name the commit actually read.
2. Write the body to a file. It opens with a note naming the coding agent, because `gh`
   posts as the person who is logged in and the review would otherwise read as theirs. The marker comes next, exactly this shape
   (full 40-character SHA, date in YYYY-MM-DD):
   ```
   > [!NOTE]
   > Written by <coding agent> (<model>) at <github user>'s request. Posted from their account, so replies here are to the agent, not to them.

   <!-- rustacean-review sha=<head sha> date=<today> -->
   ## /rustacean-review, run <YYYY-MM-DD HH:MM> UTC

   **Commit:** `<short sha>`, the PR's head when the review was read.

   <the review, in the Output format above>
   ```
   The heading, time and commit are what Olga reads on the PR to see that a review happened
   and which commit it covered. The hidden marker line only feeds the weekly table. Never
   post the marker without the visible heading and the review under it.
3. Post it with `gh api` (never `gh pr comment` or `gh pr edit`), and read it back to check
   that it starts with the note and carries the marker:
   ```bash
   gh api -X POST repos/<owner>/<repo>/issues/<n>/comments -F body=@review.md --jq .html_url
   gh api repos/<owner>/<repo>/issues/<n>/comments --jq '.[-1].body' | head -4
   ```
4. Give Olga the comment link. When a review covers a stack of PRs, post one comment per PR,
   each covering that PR's own diff and naming its own head SHA.

---

## Honesty rules

These are what make the review worth reading:

- **Separate measured from suspected.** "This allocates per iteration" is an observation.
  "This is 3× slower" is a claim requiring a benchmark. Never fabricate the second.
- **State unverified assumptions.** "If `input` can be empty, this slice panics — I couldn't
  find a caller that guarantees non-empty."
- **Never invent an API.** If you didn't look it up, don't cite it.
- **Say when code is good.** A review with an empty "Blocking" section and a specific
  "Verified clean" section is more trustworthy than one that manufactures concerns.
- **Don't cargo-cult.** "Avoid clone" is not a rule; "this clone copies 40 MB per call" is a
  finding. If you can't say why a pattern is worse here, don't flag it.
- **Push back on the premise when warranted.** If the right review comment is "this whole
  module could be twenty lines using `itertools::group_by`," say that instead of polishing
  the existing structure.

---

## Reference files

Read these when the review touches the relevant area — they hold the detailed catalogs so
this file stays scannable.

- `allocation-and-iteration.md` — before/after pairs for every clone, allocation,
  and redundant-pass pattern above, plus the borrow-restructuring recipes.
- `correctness-and-style.md` — panics, overflow, async/concurrency, `unsafe`,
  error-type design, and the API Guidelines naming table.
- `what-olga-flags.md` — Olga's own review comments on kmerseek PRs #36 to #75, grouped
  by topic, with the command to refresh the list. Read it for every kmerseek review.
- `verifying-apis.md` — exact commands and URL patterns for checking std and
  crate APIs against the pinned version, plus a stabilization-version cheat sheet.

## Before the verdict: are you sure?

Run the `are-you-sure` pass on the code under review and on any fix you wrote (added
2026-09-30). For Rust, also:
- Write or run a test for each Blocking finding that fails before the fix and passes after.
- Run `cargo test` and `cargo clippy --all-targets` on the pinned toolchain, and report both.
- For a performance claim, measure it (a benchmark or `time` on a real index), do not reason it.
