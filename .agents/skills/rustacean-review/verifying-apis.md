# Verifying APIs instead of recalling them

The single most damaging thing a review can do is confidently cite an API that doesn't exist,
was renamed, or isn't available on the project's toolchain. Rust ships every six weeks and
crates break in major versions. Treat recalled signatures as hypotheses.

## Contents

- [Establish the toolchain and versions](#establish-the-toolchain-and-versions)
- [Checking a crate API](#checking-a-crate-api)
- [Checking a std API](#checking-a-std-api)
- [MSRV gotchas](#msrv-gotchas)
- [How to write about something you couldn't verify](#how-to-write-about-something-you-couldnt-verify)

---

## Establish the toolchain and versions

Do this once, at the start of every review.

```bash
rustc --version
cargo --version
cat rust-toolchain.toml 2>/dev/null      # pins the whole team to one toolchain
grep -E '^(edition|rust-version)' Cargo.toml
```

`rust-version` in `Cargo.toml` is the MSRV. It is the ceiling on what you may suggest. If it
is absent, the project has no stated floor — ask which toolchains they support rather than
assuming the newest.

Then get the exact version of any dependency under discussion:

```bash
cargo tree -i <crate>                    # who depends on it, and at what version
grep -A1 'name = "<crate>"' Cargo.lock
```

`Cargo.toml` says `serde = "1.0"`; `Cargo.lock` says which 1.0.x is actually compiled. Read
the lock file.

---

## Checking a crate API

In descending order of reliability:

**1. Vendored source — the ground truth for the exact build.**

```bash
ls ~/.cargo/registry/src/*/            # find the crate-version directory
rg 'pub fn <name>' ~/.cargo/registry/src/*/<crate>-<version>/src/
```

If the crate isn't vendored yet, `cargo fetch` pulls it.

**2. Locally generated docs.**

```bash
cargo doc -p <crate> --no-deps
# then read target/doc/<crate>/index.html, or grep the generated HTML
```

**3. docs.rs, with the version pinned in the URL.**

```
https://docs.rs/<crate>/<exact-version>/<crate>/
https://docs.rs/<crate>/<exact-version>/<crate>/struct.<Type>.html
```

Never fetch `https://docs.rs/<crate>/latest/` when reviewing — latest may be a major version
ahead of what they compile.

**4. The crate's own changelog** for anything that looks renamed:

```
https://github.com/<org>/<repo>/blob/main/CHANGELOG.md
```

Also worth checking `Cargo.toml`'s `features` table — a method may exist but be gated behind
a feature the project hasn't enabled. `cargo tree -f '{p} {f}'` shows enabled features.

---

## Checking a std API

```bash
rustup doc --std                         # offline, matches the installed toolchain
```

Online: `https://doc.rust-lang.org/std/` for current stable, or
`https://doc.rust-lang.org/1.XX.0/std/` to check what existed at their MSRV.

Every std item's docs carry a "since 1.XX.0" stability badge — read it before suggesting the
item. Release notes at `https://github.com/rust-lang/rust/blob/master/RELEASES.md` list
stabilizations per version.

The fastest check of all is often just to write the suggestion into a scratch file and
compile it:

```bash
cargo check --all-targets
```

If you are proposing a non-trivial rewrite, do this. A suggestion that doesn't compile costs
the author more time than no suggestion at all.

---

## MSRV gotchas

Things that commonly appear in reviews and are newer than a conservative project's MSRV.
**Verify each against the current docs rather than trusting this table** — it is a reminder
that version-gating exists, not an authority.

- `let ... else` statements
- `if let` chains (`if let Some(x) = a && ...`)
- async fn in traits
- `impl Trait` in associated type position / RPITIT
- C-string literals and other edition-2024-gated syntax
- Many `Option`/`Result` combinators (`is_none_or`, `inspect`, `is_some_and`)
- `std::thread::scope`
- `OnceLock` / `LazyLock` (versus the older `once_cell` crate)
- Const generics beyond the minimal subset
- Newer `slice` methods (`array_chunks`, `split_at_checked`, and friends)

When the useful suggestion is above their MSRV, say both parts: "`LazyLock` is the clean
version of this, stable since 1.XX — above your stated MSRV of 1.YY, so either `once_cell`
or bumping the floor."

---

## How to write about something you couldn't verify

Do not stay silent about a real concern just because you couldn't confirm the fix. State the
observation as fact and the fix as a lead:

> This re-hashes `key` on both the lookup and the insert. The `entry` API collapses that
> into one hash — worth checking whether the borrow works out here given `key` is used again
> below.

Never do this:

> Use `HashMap::get_or_insert_with_ref` here.

If you didn't look it up, don't name it. A reviewer who invents one API loses the author's
trust in all the rest of the findings.
