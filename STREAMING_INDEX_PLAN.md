# Streaming index: removing the whole-corpus memory ceiling

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

Peak RSS is ~500 bytes per residue and scales linearly. UniRef50 is roughly
70 M sequences / ~18 G residues, about 900x that subset, which extrapolates to
**~8 TB of RAM per encoding** — ~65x the available 128 GB. Reducing the
alphabet does not help: hp k=24 is only 16% cheaper than protein k=9, because
the cost is dominated by per-sequence storage rather than by k-mer diversity.

## The five O(corpus) structures held in RAM

All verified by reading `src/rust/index.rs` on `main` (c1d8135).

1. **`signatures: DashMap<String, ProteinSketch>`** (index.rs:91)
   Every sketch is retained for the whole run. `process_fasta` already streams
   the FASTA in batches, but `process_batch_parallel` inserts into this map and
   nothing is ever evicted. Per sketch the dominant field is
   `kmer_positions: HashMap<u64, Vec<usize>>`.

2. **`rebuild_combined_minhash()`** (index.rs:363)
   Collects *every* min from *every* signature into one `Vec<u64>` before
   `sort_unstable` + `dedup`. This is pre-dedup, so it holds one u64 per k-mer
   *instance*, not per unique k-mer — ~18 G entries (~144 GB) for UniRef50, in a
   single contiguous allocation.

3. **`save_state()` → `signature_data: Vec<ProteinSketchStore>`** (index.rs:396)
   Converts all signatures into a second full in-memory copy before chunking.

4. **`ProteomeIndexMetadata.combined_mins: Vec<u64>`**
   Serialized to one value. Measured at **859 MB for only 250 K sequences**;
   linear in unique k-mers.

5. **`save_inverted_index()`**
   Builds `inverted_index: HashMap<u64, Vec<u32>>` and
   `kmer_frequencies: HashMap<u64, usize>` fully in RAM, then bincodes the whole
   `SearchCache` into **one contiguous `Vec<u8>` for a single `db.put`**. This
   alone caps an index at what one allocation can hold.

## Direction

Ordered by expected benefit per unit of risk.

1. **Evict signatures after each batch.** Write each `ProteinSketchStore` to its
   `sig_{md5}` key as its batch completes and drop it from the DashMap. This is
   the single biggest win and is mostly independent of the rest. Requires
   deciding what search still needs resident (currently `sig_cache` handles
   on-demand loads, so possibly nothing).

2. **Accumulate the inverted index in RocksDB, not RAM.** One key per k-mer with
   a merge operator, or fixed hash-range shards flushed periodically. Removes
   both hotspot 5 and the single-blob `SearchCache` write.

3. **Derive the combined minhash from the inverted index keys.** The set of
   unique k-mers *is* the inverted index key set, so hotspot 2's giant pre-dedup
   Vec is redundant work. Feed hashes in sorted order (already the cheap path
   per the existing `rebuild_combined_minhash` comment).

4. **Drop `signature_data`.** Serialize each signature straight into its chunk
   rather than materializing the full Vec first.

5. **Store `combined_mins` chunked** (or stop storing it, if 3 makes it
   derivable at open time).

## Validating

- Peak RSS must become roughly flat in corpus size rather than linear. Measure
  with `/usr/bin/time -l` across increasing SwissProt subsets (50 K / 250 K /
  573 K sequences) and compare the slope against the table above.
- Indexes built before and after must be equivalent: existing
  `test_index_equivalence` and `test_manual_vs_auto_index_equivalence` cover
  this, plus search results on the bcl2/ced9 fixtures must be unchanged.
- `SCHEMA_VERSION` needs a bump if the on-disk layout changes (it will for
  items 2 and 5).
