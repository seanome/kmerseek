# Allocation, cloning, and iteration patterns

Before/after pairs for the findings in Phase 2 and Phase 3. Every "after" must be checked
against the crate's MSRV before you suggest it.

## Contents

- [Clone triage](#clone-triage)
- [Restructuring instead of cloning](#restructuring-instead-of-cloning)
- [Allocation in loops](#allocation-in-loops)
- [String and buffer building](#string-and-buffer-building)
- [Parameter and return types](#parameter-and-return-types)
- [Redundant passes](#redundant-passes)
- [Map and set access patterns](#map-and-set-access-patterns)
- [Quadratic patterns](#quadratic-patterns)
- [Parallelism](#parallelism)
- [What not to flag](#what-not-to-flag)

---

## Clone triage

### Clone on a `Copy` type

```rust
let n = count.clone();   // count: usize
let n = count;           // Copy is a memcpy; clone() adds nothing but confusion
```

Clippy: `clone_on_copy`. Free fix, low importance — group these into one line of the review
rather than one finding each.

### Cloning a whole struct to read one field

```rust
let name = record.clone().name;      // deep-copies every field, drops all but one
let name = &record.name;             // or record.name.clone() if ownership is needed
```

### Cloning a key for lookup

```rust
if map.contains_key(&key.clone()) { ... }
if map.contains_key(&key) { ... }

// String keys accept &str via Borrow — no allocation needed
let v = map.get(&owned_string);      // fine
let v = map.get(name_str);           // also fine, map: HashMap<String, _>
```

### `Arc`/`Rc` clone — this one is fine

```rust
let handle = Arc::clone(&shared);    // prefer this spelling
let handle = shared.clone();         // works, but reads like a deep copy
```

Refcount bumps are the intended use. Do not flag them as waste. The only note worth making
is the spelling: `Arc::clone(&x)` tells the reader at a glance that nothing was duplicated.

### Cloning to satisfy a `'static` bound

```rust
thread::spawn(move || process(data.clone()));
```

Often correct — the thread genuinely needs ownership. Check whether `Arc` or a scoped thread
(`std::thread::scope`, stable since 1.63) removes the need to duplicate the data at all.

---

## Restructuring instead of cloning

When a clone exists only to end a borrow, the fix is structural.

### Split the borrow

```rust
// before: cannot borrow `self` mutably while `self.config` is borrowed
let cfg = self.config.clone();
self.update(&cfg);

// after: destructure so the two fields are borrowed independently
let Self { config, state, .. } = self;
state.update(config);
```

### `std::mem::take` / `std::mem::replace`

```rust
// before: clone to move data out of &mut self
let items = self.queue.clone();
self.queue.clear();

// after: one move, no copy
let items = std::mem::take(&mut self.queue);
```

### Index instead of holding a reference

```rust
// before: clone because the element borrow conflicts with the mutation
let item = list[i].clone();
list.push(transform(&item));

// after: compute first, then mutate
let new = transform(&list[i]);
list.push(new);
```

### Reorder so the borrow ends first

Non-lexical lifetimes mean a borrow ends at its last use, not at the end of the block. Often
moving one statement up removes the conflict entirely.

---

## Allocation in loops

The highest-value category in data-processing code.

### Hoist and reuse a scratch buffer

```rust
// before: one Vec per record
for record in records {
    let mut buf = Vec::new();
    encode_into(record, &mut buf);
    sink.write_all(&buf)?;
}

// after: one Vec total
let mut buf = Vec::new();
for record in records {
    buf.clear();                 // keeps the capacity, drops the contents
    encode_into(record, &mut buf);
    sink.write_all(&buf)?;
}
```

Note `clear()` retains capacity — that is the whole point. If the buffer can grow unbounded
across iterations, add a `shrink_to` guard rather than reverting to per-iteration allocation.

### Preallocate when the length is known

```rust
let mut out = Vec::new();                    // log2(n) reallocs + copies
let mut out = Vec::with_capacity(items.len());
```

`Iterator::collect` already does this when the iterator has an exact size hint, so
`items.iter().map(f).collect()` needs no help. The manual push loop does.

### `extend` beats a push loop

```rust
for x in other { v.push(x); }
v.extend(other);                 // uses the size hint to reserve once
```

---

## String and buffer building

```rust
// before: allocates a String per iteration, then another for the join
let parts: Vec<String> = ids.iter().map(|i| format!("{i}")).collect();
let line = parts.join(",");

// after: one buffer
use std::fmt::Write;
let mut line = String::with_capacity(ids.len() * 8);
for (n, id) in ids.iter().enumerate() {
    if n > 0 { line.push(','); }
    write!(line, "{id}").unwrap();   // writing to a String is infallible
}
```

Also: `s.push_str(x)` over `s = s + x`; `s.push('c')` over `s.push_str("c")`.

### `Cow` for the mostly-borrowed case

```rust
fn normalize(s: &str) -> Cow<'_, str> {
    if s.chars().all(|c| c.is_ascii_lowercase()) {
        Cow::Borrowed(s)               // common case: no allocation
    } else {
        Cow::Owned(s.to_lowercase())
    }
}
```

Worth suggesting only when the borrowed case genuinely dominates. Otherwise it adds a match
at every use site for nothing.

---

## Parameter and return types

| Instead of | Take | Why |
|---|---|---|
| `&Vec<T>` | `&[T]` | Accepts arrays, slices, `SmallVec`, and `Vec` |
| `&String` | `&str` | Accepts literals and slices without allocation |
| `&PathBuf` | `&Path` or `impl AsRef<Path>` | Callers stop building a `PathBuf` just to call you |
| `String` (when only read) | `&str` | Stops forcing the caller to clone |
| returning `&Vec<T>` | returning `&[T]` | Doesn't leak the container choice into the API |
| `Vec<T>` param you only iterate | `impl IntoIterator<Item = T>` | Caller need not collect first |

Counterpoint: take `String` by value when you *will* store it — `impl Into<String>` is the
polite version, since callers holding a `String` pay nothing and callers holding a `&str`
allocate exactly once.

---

## Redundant passes

### Multiple statements over the same collection

```rust
// before: three traversals
let total = v.iter().filter(|x| x.ok).count();
let sum: u64 = v.iter().filter(|x| x.ok).map(|x| x.len).sum();
let max = v.iter().filter(|x| x.ok).map(|x| x.len).max();

// after: one
let (total, sum, max) = v.iter().filter(|x| x.ok).fold(
    (0usize, 0u64, None::<u64>),
    |(n, s, m), x| (n + 1, s + x.len, m.max(Some(x.len))),
);
```

Judgment call: if the fold becomes unreadable, three clear passes over a small collection is
the better code. Flag it only when the collection is large or the traversal is expensive.

### Mid-chain `collect`

```rust
let names: Vec<String> = users.iter().map(|u| u.name.clone()).collect();
let count = names.iter().filter(|n| n.starts_with('a')).count();

let count = users.iter().filter(|u| u.name.starts_with('a')).count();
```

Clippy: `needless_collect`. Note the clone disappeared too — that is typical, the collect
was forcing ownership that was never needed.

### `collect` then immediately consume

```rust
let v: Vec<_> = iter.collect();
for x in v { ... }

for x in iter { ... }
```

---

## Map and set access patterns

### Double hashing

```rust
// before: hashes `key` twice on the insert path
if !map.contains_key(&key) {
    map.insert(key.clone(), default());
}
let v = map.get_mut(&key).unwrap();     // and a third time

// after
let v = map.entry(key).or_insert_with(default);
```

`or_insert_with` over `or_insert` when the default is expensive to construct — `or_insert`
evaluates its argument unconditionally.

### Accumulating into a map

```rust
*counts.entry(k).or_insert(0) += 1;
counts.entry(k).or_default().push(v);
```

### Choosing the container

- `HashMap` — the default. Consider `rustc_hash::FxHashMap` or `ahash` when keys are small
  and untrusted input isn't a concern; SipHash is DoS-resistant but not fast.
- `BTreeMap` — when you need ordered iteration or range queries.
- Sorted `Vec` + `binary_search` — best cache behavior for build-once/query-many, and often
  faster than a hash map for small n.
- `HashSet` vs sort+dedup — `HashSet` when order is irrelevant and n is large; sort+dedup
  when you need the sorted output anyway or want to avoid hashing cost.

---

## Quadratic patterns

### Membership check inside a loop

```rust
// before: O(n·m)
for x in &a {
    if b.contains(x) { ... }        // b: &Vec<_>, linear scan each time
}

// after: O(n + m)
let b: HashSet<_> = b.iter().collect();
for x in &a {
    if b.contains(x) { ... }
}
```

### Removing inside a loop

```rust
// before: each remove shifts the tail — O(n²)
let mut i = 0;
while i < v.len() {
    if !keep(&v[i]) { v.remove(i); } else { i += 1; }
}

// after: one pass
v.retain(keep);
```

`swap_remove` is O(1) if element order doesn't matter.

### Repeated `insert(0, ..)` or front removal

`Vec` is not a deque. Use `VecDeque` when both ends see traffic.

### String concatenation in a loop without capacity

Each `push_str` past capacity reallocates and copies. `String::with_capacity` up front, or
`write!` into a preallocated buffer.

---

## Parallelism

Before approving or suggesting `rayon`:

- **Work per item must exceed the scheduling overhead.** `par_iter()` over 100 integers doing
  an addition is slower than the serial loop.
- **`par_bridge` is not free** — it serializes through a shared iterator. Prefer a source
  that implements `IntoParallelIterator` directly.
- **Contended shared state kills the win.** A `Mutex<Vec<_>>` written by every task means the
  work is effectively serial plus lock overhead. Use `fold`/`reduce` or `collect` into a
  per-thread accumulator.
- **`chunks(n).par_bridge()`** amortizes overhead when per-item work is small but total work
  is large.
- Any parallelism claim needs a benchmark. Say so.

---

## What not to flag

Reviews lose credibility fast when they include these:

- Clones in `main`, `build.rs`, tests, benchmarks, CLI argument parsing, or any path that
  runs once. The cost is unmeasurable and the clarity is real.
- Iterator adapter chains — they fuse; there is no extra pass.
- `Arc::clone` on a shared handle.
- `to_string()` in an error path that already allocates a `Box<dyn Error>`.
- `collect()` where the result is genuinely stored or returned.
- Micro-optimizations in code that is not hot, unless the idiomatic version is also shorter.
- Anything you would have to guess about. Ask instead.
