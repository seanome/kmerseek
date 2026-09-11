# Streaming index: removing the whole-corpus memory ceiling

Status: implemented in this branch. The sections below keep the measurements
and reasoning that motivated it; "What was done" at the end says how each
structure was removed and what the result measures.

## Why

Indexing currently holds the entire corpus in RAM. That caps the largest
indexable proteome well below UniRef50.

Measured on this machine (128 GB RAM), release build, three encodings over an
identical 50,000-sequence SwissProt subset (19.7 M residues):

| Encoding      | Unique k-mers | Peak RSS | Disk   | Time   |
|---------------|---------------|----------|--------|--------|
| protein k=9   | 11.15 M       | 9.8 GB   | 1.2 GB | 14.4 s |
| dayhoff k=17  | 10.54 M       | 9.1 GB   | 1.1 GB | 13.5 s |
| hp k=24       | 7.32 M        | 8.2 GB   | 1.0 GB | 11.8 s |

Peak RSS is ~500 bytes per residue and scales linearly. A later measurement on
a UniRef50 sample at k=10 protein20 gave 683 B/residue (Swiss-Prot control:
634), so the coefficient is data-independent to within 8%. UniRef50 release
2026_03 is 38.8 M sequences / 12.18 G residues, about 600x that subset, which
extrapolates to ~8.3 TB of RAM per encoding, 65x the available 128 GB.

Reducing the alphabet does not help much: hp k=24 is only 16% cheaper than
protein k=9. The reason is that the cost is per k-mer *instance*, not per unique
k-mer, and every alphabet produces one instance per residue position. The
`--scaled` measurements in PR #50 confirm this: RSS per residue fits
~21 + 662/scaled, so about 97% of indexing memory is per-k-mer and only ~3% is
per-sequence overhead. FracMinHash sampling (`--scaled`) therefore lowers the
constant almost linearly, but leaves the O(corpus) scaling in place. This plan
is about removing the scaling; `--scaled` is the complementary lever for the
constant.

## The five O(corpus) structures held in RAM

All verified by reading `src/rust/index.rs` on `main` (5f5b0b8).

1. **`signatures: DashMap<String, ProteinSketch>`** (index.rs:140)
   Every sketch is retained for the whole run. `process_fasta` already streams
   the FASTA in batches, but `process_batch_parallel` inserts into this map and
   nothing is ever evicted. Per sketch the dominant field is
   `kmer_positions: HashMap<u64, Vec<usize>>`.

2. **`rebuild_combined_minhash()`** (index.rs:884)
   Collects *every* min from *every* signature into one `Vec<u64>` before
   `sort_unstable` + `dedup`. This is pre-dedup, so it holds one u64 per k-mer
   *instance*, not per unique k-mer: ~12 G entries (~97 GB) for UniRef50, in a
   single contiguous allocation.

3. **`save_state_with_kmer_stats()` → `signature_data: Vec<ProteinSketchStore>`** (index.rs:997)
   Converts all signatures into a second full in-memory copy before chunking.

4. **`ProteomeIndexMetadata.combined_mins` and `combined_abunds`** (index.rs:71)
   Both `Vec<u64>`, serialized together into the single `index_metadata` value,
   which measures exactly 16·U + 60 bytes (U = unique k-mers). The abundances
   are there because `rebuild_combined_minhash` passes `track_abundance = true`.
   Measured at 859 MB for only 250 K sequences; Swiss-Prot at k=10 already sits
   at 43% of RocksDB's 4 GiB single-value limit, which UniRef50 exceeds by 42x.
   Nothing on the search path reads either field: `open_for_search` (index.rs:1373)
   builds an empty combined minhash and says so, and the only readers outside
   `save_state` are tests. Search startup still pays to deserialize the whole
   value.

5. **`save_inverted_index()`** (index.rs:442)
   Builds `inverted_index: HashMap<u64, Vec<u32>>` and
   `kmer_frequencies: HashMap<u64, usize>` fully in RAM, then bincodes the whole
   `SearchCache` (40·T + 32·U + 4·P bytes for T targets, U unique k-mers, P
   postings) into one contiguous `Vec<u8>` for a single `db.put`. PR #48 chunks
   that write, which lifts the 4 GiB value limit; the in-RAM build is untouched.

## Direction

Ordered by expected benefit per unit of risk.

1. **Evict signatures after each batch.** Write each `ProteinSketchStore` to its
   `sig_{md5}` key as its batch completes and drop it from the DashMap. This is
   the single biggest win and is mostly independent of the rest. Requires
   deciding what search still needs resident (currently `sig_cache` handles
   on-demand loads, so possibly nothing).

2. **Accumulate the inverted index in RocksDB, not RAM.** One key per k-mer with
   a merge operator, or fixed hash-range shards flushed periodically. Removes
   hotspot 5's in-RAM build and supersedes the chunked write from PR #48.

3. **Derive the combined minhash from the inverted index keys.** The set of
   unique k-mers *is* the inverted index key set, so hotspot 2's giant pre-dedup
   Vec is redundant work. Feed hashes in sorted order (already the cheap path
   per the existing `rebuild_combined_minhash` comment).

4. **Drop `signature_data`.** Serialize each signature straight into its chunk
   rather than materializing the full Vec first.

5. **Stop storing `combined_mins` and `combined_abunds`.** Nothing on the
   search path reads them (hotspot 4), and item 3 makes the mins derivable at
   open time for the tests that do. This also removes the next hard format
   limit after PR #48, the 4 GiB `index_metadata` value.

## Validating

- Peak RSS must become roughly flat in corpus size rather than linear. Measure
  with `/usr/bin/time -l` across increasing SwissProt subsets (50 K / 250 K /
  573 K sequences) and compare the slope against the table above.
- Indexes built before and after must be equivalent: existing
  `test_index_equivalence` and `test_manual_vs_auto_index_equivalence` cover
  this, plus search results on the bcl2/ced9 fixtures must be unchanged.
- `SCHEMA_VERSION` needs a bump if the on-disk layout changes (it will for
  items 2 and 5).
- Peak RSS at `--scaled 1` must stay at or below the measured 683 B/residue
  at every step, so a partial implementation does not regress the constant
  while attacking the slope.

## What was done

Schema version 3. `src/rust/index.rs` no longer holds anything O(corpus):

1. **Signatures stream to disk.** `ingest` writes each sketch to its `sig_{md5}`
   key in one `WriteBatch` per batch and drops it. The in-memory map is now
   only for `store_signatures` and `load`, whose callers want every sketch in
   hand. The dedup the map used to do by key is kept with a `HashSet<u64>`
   (8 bytes per sequence); the first occurrence wins, so which name a repeated
   sequence carries is now fixed by FASTA order rather than by which parallel
   insert landed last.

2. **The inverted index is built by an external sort.** `(hash, target)` pairs
   go into a buffer of 16 M entries (256 MB). When it fills it is bucketed by
   `hash & 1023` and written as one run per shard (`ii_run_{r}_{s}`).
   `finalize` reads each shard's runs, sorts by `(hash, target)`, groups into
   posting lists, writes `ii_shard_{s}`, and deletes the runs. Peak memory is
   one shard, about 1/1024 of the postings. Shards use the low bits because
   FracMinHash keeps only hashes below `max_hash / scaled`, which would leave
   the high bits nearly constant.

3. **The combined minhash is gone.** `unique_kmer_count()` is the number of
   inverted index keys, counted while the shards are written.

4. **`signature_data` and the `signatures_chunk_{n}` keys are gone.** `load`
   reads the `sig_{md5}` keys with a prefix iterator.

5. **`combined_mins` and `combined_abunds` are gone from the metadata**, so the
   next 4 GiB single-value limit is gone with them. `kmer_frequencies` is no
   longer stored either: it is the length of each posting list, derived when
   the search cache is loaded.

The k-mer frequency spectrum and the most/least common examples are gathered in
the same pass that writes the shards, with two bounded heaps (`SmallestN`), so
`--kmer-stats-out` needs no frequency map. The examples' text is recovered by
reading the first target of each from disk. `--stats-only` builds in a scratch
directory that is removed on exit, since the sort needs the disk.

The writer no longer opens RocksDB with mmap. Every SST page touched during
compaction was otherwise mapped into the process and counted in RSS, which left
a slope of ~70 bytes per residue, the on-disk size, after the in-memory
structures were removed.

Indexes written by schema 2 still open: `load_search_cache` reads their single
`search_cache` value and `load` their signature chunks. PR #48's chunked cache
write is superseded.

### Measured

Swiss-Prot 2026_03 prefixes and the full release, k=10 protein20, release
builds, `/usr/bin/time -l`, same machine as the table above.

| Sequences | Residues | Old peak RSS | New peak RSS | Old time | New time | Old disk | New disk |
|-----------|----------|-------------:|-------------:|---------:|---------:|---------:|---------:|
| 20,000    | 8.3 M    | 4.88 GB      | 1.01 GB      | 7.4 s    | 2.1 s    | 532 MB   | 255 MB   |
| 50,000    | 19.6 M   | 10.35 GB     | 1.82 GB      | 17.7 s   | 4.7 s    | 1.22 GB  | 587 MB   |
| 100,000   | 38.7 M   | 19.35 GB     | 2.21 GB      | 29.3 s   | 10.1 s   | 2.42 GB  | 1.17 GB  |
| 575,748   | ~200 M   | not run (~100 GB by the old slope) | 2.88 GB | | 108 s | | 6.14 GB |

The old code is 500 B/residue. The new code's growth from 100 K sequences to
the full release is 0.67 GB over ~160 M residues, about 4 B/residue, and most
of that is the posting buffer and RocksDB memtables filling up to their fixed
sizes. Disk halves because signatures are stored once instead of twice.

Searching the 25 BCL2 query sequences against the 50 K index gives the same
110 hits and 458 rows from the old index, the new index, and the old index read
by the new binary; the only numeric differences are at 1e-15 (summation order).
