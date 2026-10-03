# Correctness and style

Detail for Phase 1 and Phase 4. Check every suggested API against the crate's MSRV before
recommending it — see `verifying-apis.md`.

## Contents

- [Panics](#panics)
- [Arithmetic and casts](#arithmetic-and-casts)
- [Error design](#error-design)
- [Unsafe](#unsafe)
- [Trait contracts](#trait-contracts)
- [Concurrency and async](#concurrency-and-async)
- [Naming conventions](#naming-conventions)
- [API shape](#api-shape)
- [Module and visibility structure](#module-and-visibility-structure)
- [Tests](#tests)

---

## Panics

Every panic is an undocumented assertion. In a binary that may be fine; in a library it is a
contract the caller cannot see.

| Pattern | Question to ask |
|---|---|
| `.unwrap()` | Is the invariant guaranteed by a caller, a type, or nothing? |
| `.expect("...")` | Does the message say *why it can't happen*, not *what failed*? Good: `"k is validated non-zero at construction"`. Bad: `"unwrap failed"`. |
| `v[i]` | Is `i` derived from `v`'s own length, or from user input? `.get(i)` returns `Option`. |
| `&s[a..b]` | On a `str`, this panics on non-char-boundary — a real bug with UTF-8 input. |
| `a / b`, `a % b` | Can `b` be zero? |
| `.unwrap()` on `Mutex::lock` | Conventional and usually fine — poisoning means another thread already panicked. Worth a comment, not a rewrite. |

Better shapes: `let Some(x) = opt else { return Err(...) };`, `ok_or_else`, `?`. If a panic
is genuinely the right behavior, document it under `# Panics` and prefer `expect` with a
reason over bare `unwrap`.

Also flag: `unwrap()` in a `Drop` impl (double panic aborts), and `assert!` used for input
validation on a public API where an error would serve the caller better.

---

## Arithmetic and casts

Debug builds check overflow; **release builds wrap silently by default**. Code tested only
in debug can be wrong in production.

```rust
let n = a + b;                       // wraps in release
let n = a.checked_add(b).ok_or(Error::Overflow)?;
let n = a.saturating_add(b);
let n = a.wrapping_add(b);           // fine — but say so, so readers know it's deliberate
```

Casts:

```rust
let small = big as u32;              // silently truncates
let small = u32::try_from(big)?;     // fails loudly

let n = x as usize;                  // on a negative i64 this becomes enormous
let f = big_u64 as f64;              // loses precision above 2^53
let i = f as i32;                    // saturates since 1.45 — was UB before; check MSRV
```

Flag any `as` between integer types of different width or signedness unless a comment
explains why the truncation is intended.

---

## Error design

**Libraries**: a concrete error enum callers can match on.

```rust
#[derive(Debug, thiserror::Error)]
pub enum SketchError {
    #[error("k must be greater than zero")]
    InvalidK,
    #[error("failed to read {path}")]
    Io { path: PathBuf, #[source] source: std::io::Error },
}
```

Keep the `#[source]` chain intact — that is what gives the caller a useful backtrace.

**Binaries**: `anyhow` with `.context(...)` at each layer.

Anti-patterns:

- `Box<dyn Error>` or `String` in a public library signature — the caller can only print it.
- `map_err(|_| MyError::Generic)` — discards the cause.
- `.ok()` or `let _ = fallible()` with no comment — silent failure.
- An error enum with a single `Other(String)` variant that everything funnels into.
- `unwrap_or_default()` where the default is indistinguishable from a real value.

---

## Unsafe

Every `unsafe` block needs a `// SAFETY:` comment naming the specific invariant being upheld
and why it holds here. Not "this is safe" — *which* precondition, and *what* guarantees it.

```rust
// SAFETY: `idx < self.len` is checked on the line above, and `self.ptr` is valid for
// `self.len` elements for the lifetime of `self` (established in `Self::new`).
let value = unsafe { *self.ptr.add(idx) };
```

Review checklist:

- Could safe code do this? Measure before accepting `unsafe` for speed.
- `from_utf8_unchecked` — is the input *guaranteed* UTF-8, or just usually?
- `get_unchecked` — is the bound checked immediately above, or several frames away?
- Raw pointer arithmetic — provenance, alignment, and null.
- `transmute` — almost always replaceable with a safe cast, `bytemuck`, or `zerocopy`.
- Does the surrounding safe API make it impossible for a user to trigger UB? An `unsafe`
  block in a `pub fn` with no `unsafe` marker means the *function* must uphold the invariant.
- Recommend `cargo +nightly miri test` for anything touching raw pointers.

---

## Trait contracts

- `Hash` and `Eq` must agree: `a == b` implies `hash(a) == hash(b)`. Deriving both is safe;
  hand-implementing one and deriving the other is a classic bug.
- `Ord` must be consistent with `PartialOrd` and total. A hand-rolled `Ord` that returns
  `Equal` for incomparable values silently corrupts `sort` and `BTreeMap`.
- Floats are not `Ord` for a reason. `f64` as a `HashMap` key or in a `BTreeSet` is a smell —
  `ordered-float` or a fixed-point integer is usually the right answer.
- `PartialEq` that isn't symmetric or transitive.
- `Default` that produces an invalid instance of a type with invariants.
- `Clone` hand-implemented without matching `clone_from`, when `clone_from` could reuse
  allocations.
- Code that depends on `HashMap`/`HashSet` iteration order — it is unspecified and varies
  per run.

---

## Concurrency and async

- **`MutexGuard` held across `.await`** — `std::sync::MutexGuard` is not `Send`, so this
  either fails to compile or (with `tokio::sync::Mutex`) silently serializes the executor.
  Restructure to drop the guard before awaiting.
- **Blocking I/O or CPU work in an async fn** — starves the runtime. `spawn_blocking`.
- **`Rc`/`RefCell` in code that must be `Send`** — should be `Arc`/`Mutex`.
- **Lock ordering** — two locks acquired in different orders in different functions is a
  deadlock waiting for load.
- **`Ordering::Relaxed`** on an atomic used to publish other data. `Relaxed` orders nothing
  but the atomic itself; `Acquire`/`Release` is what makes the adjacent writes visible.
- **Unbounded channels** as backpressure-free queues — memory grows without limit.
- **Detached tasks** whose `JoinHandle` is dropped — errors vanish silently.
- **`std::sync::mpsc` vs `crossbeam`/`flume`** — worth mentioning only if the limitation bites.

---

## Naming conventions

From the Rust API Guidelines. These names carry cost information; getting them backwards
misleads every caller.

| Prefix | Cost | Receiver | Example |
|---|---|---|---|
| `as_` | Free | borrowed → borrowed | `str::as_bytes` |
| `to_` | Expensive | borrowed → owned | `str::to_string` |
| `into_` | Variable | owned → owned, consumes `self` | `String::into_bytes` |

Other conventions worth enforcing:

- Getters have no `get_` prefix: `fn name(&self) -> &str`, not `fn get_name`.
- `iter` / `iter_mut` / `into_iter` for the three iterator flavors.
- `is_` / `has_` for predicates returning `bool`.
- Types are `UpperCamelCase`, functions and variables `snake_case`, constants `SCREAMING_SNAKE`.
- Acronyms are one word in camel case: `HttpClient`, `Uuid` — not `HTTPClient`, `UUID`.
- Don't stutter: `sketch::SketchBuilder` reads as `sketch::Builder` at the call site.

---

## API shape

- **Take generic, return concrete.** `impl AsRef<Path>`, `&str`, `&[T]`, `impl IntoIterator`
  in argument position; a named concrete type coming out.
- **Newtypes over bare primitives.** Two `usize` parameters in a row is an invitation to swap
  them at the call site. `struct KmerSize(u8)` makes that a compile error.
- **Builders** for constructors with more than ~4 arguments or many optional ones.
- **`#[non_exhaustive]`** on public enums and structs you expect to extend, so adding a
  variant isn't a breaking change.
- **`#[must_use]`** on constructors, pure transforms, and anything where discarding the
  result is certainly a bug.
- **Derive liberally**: `Debug` on every public type (its absence is a real papercut),
  plus `Clone`, `PartialEq`, `Eq`, `Hash`, `Copy`, `Default` where they genuinely fit.
  `PartialOrd`/`Ord` only when a total order actually means something for the type.
- **Trait bounds on `where` clauses** rather than inline, once there is more than one.
- **`impl Trait` return** hides the concrete type — good for iterators, but it also prevents
  callers from naming it. Consider a named type for anything long-lived in the API.

---

## Module and visibility structure

- Default to `pub(crate)`. `pub` is a semver promise.
- Prefer `mod foo;` + `foo.rs` over `foo/mod.rs` in edition 2018+.
- Re-export the public surface from the crate root so users write `crate::Thing`, not
  `crate::internal::detail::Thing`.
- Group `use` statements: std, external crates, then local — `rustfmt` with
  `group_imports = "StdExternalCrate"` (nightly-only option) does this automatically.
- A `utils.rs` or `helpers.rs` growing past a few hundred lines usually means a missing
  abstraction.

---

## Tests

- Unit tests in `#[cfg(test)] mod tests` next to the code; integration tests in `tests/`
  exercising only the public API. If the integration tests need private items, the API is
  probably wrong.
- Property tests (`proptest`, `quickcheck`) for parsers, encoders, and anything with a
  round-trip law — `decode(encode(x)) == x` catches more than a dozen hand-written cases.
- Tests asserting only that a function doesn't panic assert almost nothing.
- Benchmarks belong in `benches/` with `criterion`, not in an ad-hoc `main`. Any performance
  claim in the review should point to one.
- `#[should_panic]` without `expected = "..."` passes on the wrong panic.
