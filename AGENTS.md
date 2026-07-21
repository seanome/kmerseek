# kmerseek Project Instructions

## Project Overview
Rust protein k-mer search tool using HP/Dayhoff/protein encodings and MinHash sketches.
- Key files: `src/rust/index.rs`, `src/rust/search.rs`, `src/rust/main.rs`, `Cargo.toml`
- Branch: `olgabot/rust-search`

## Architecture
- `ProteomeIndex`: RocksDB-backed index; stores signatures in chunks + individual `sig_{md5}` entries + `search_cache`
- `ProteinSearcher`: wraps ProteomeIndex + SearchStats + target_list + inverted_index for search
- `search_one()`: streaming search path used in main.rs (one query at a time)
- `ProteinSketch`: full signature with kmer_infos needed for `find_matched_regions()`

## Performance Optimizations Implemented
1. **Inverted index** (`HashMap<u64, Vec<u32>>`): reduces O(Q×T) to O(Q×candidates)
2. **Search cache** saved at index time: `target_list`, `inverted_index`, `kmer_frequencies`
3. **Individual signature storage** (`sig_{md5}` RocksDB keys): on-demand loading at search time
4. **Fast search startup** via `ProteomeIndex::open_for_search()`: reads only metadata, no signatures loaded
5. **`sig_cache: DashMap`** in `ProteinSearcher`: lazy cache avoids repeated RocksDB reads for hot targets
6. **Batch parallel query processing** in main.rs: `BATCH_SIZE=500` queries processed via `par_iter()` at the outer loop
7. **`search()` and `search_all_vs_all()`** updated to use inverted index via `search_one()` instead of O(N²) exhaustive search
8. **PreparedQuery.mins** used directly in `compare()` instead of recomputing
9. **`calculate_similarity_from_precomputed()`**: avoids double intersection computation
10. **Sequential iter** (not par_iter) for small collections (intersection stats)

## Key API
- `ProteomeIndex::open_for_search(path)`: minimal open (no signature loading)
- `ProteomeIndex::load_search_cache()`: returns `(target_list, inverted_index, kmer_frequencies)`
- `ProteomeIndex::get_signature_by_md5(md5)`: on-demand signature loading from RocksDB
- `ProteomeIndex::save_inverted_index()`: private; called by `save_state()`
- `ProteinSearcher::load(path)`: fast path (uses cache) or slow path (loads all signatures)

## RocksDB Notes
- No column families currently (flat key namespace)
- Individual signature keys: `sig_{md5}` → serialized `ProteinSketchStore`
- Metadata: `index_metadata` key → `ProteomeIndexMetadata`
- Search cache: `search_cache` key → `SearchCache` (target_list + inverted_index + kmer_frequencies)
- Chunk format (backward compat): `signatures_chunk_{n}` keys

## Testing
- All tests: `cargo test --no-default-features --lib -- --test-threads=2`
- Single test: `cargo test --no-default-features --lib search::tests::test_search_database_bcl2_ced9`
- With default threads, occasional SIGABRT from RocksDB lock contention in parallel tests (pre-existing, not a bug in our code)
- `test_cli_search_bcl2_ced9` now checks for "Total matches" (was "Found") in stderr