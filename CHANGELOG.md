# Changelog

## 0.4.0

- Added a simple visualization of target hits against a query sequence, showing both the original and encoded alphabets ([#28](https://github.com/seanome/kmerseek/pull/28))
- Added removal of hits with uncorrected p-value > 0.05, since those are removed anyway after BH correction ([#27](https://github.com/seanome/kmerseek/pull/27))
- Added filtering by at least 2 shared hashes per hit to remove spurious hits ([#27](https://github.com/seanome/kmerseek/pull/27))
- Added a table of HP encodings to the README ([#33](https://github.com/seanome/kmerseek/pull/33))
- Removed the `--scaled` option since we always use `--scaled 1` anyway ([#30](https://github.com/seanome/kmerseek/pull/30))
- Sorted top hits in the visualization by Benjamini-Hochberg corrected p-value (q-value) instead of containment ([#32](https://github.com/seanome/kmerseek/pull/32))
- Replaced invented placeholder motifs with real CED9/BCL2 residues and their actual Thomas-Dill HP encoding in `visualize_hits` tests ([#35](https://github.com/seanome/kmerseek/pull/35))
- Logged the k-mer frequency histogram and top/bottom 10 k-mers during indexing, with `--kmer-stats-out` to save the spectrum as CSV ([#36](https://github.com/seanome/kmerseek/pull/36))
- Added `--remove-low-complexity` to drop homopolymer k-mers (e.g. poly-glutamate tracts, all-hydrophobic HP windows) at index time; the setting is stored in the index and applied automatically at search time ([#37](https://github.com/seanome/kmerseek/pull/37))
- Added region-scoped scoring (`region_poisson_score`, `region_tail_probability`, `region_enrichment`, and related columns) so a matched sub-region is scored on its own instead of being diluted by the whole protein; added `--min-region-score`, deprecated `--max-pvalue` in favor of `--max-query-pvalue`, and renamed several result columns (`poisson_pvalue` → `query_poisson_pvalue`, `query_start`/`query_end`/`query_subseq` → `region_start`/`region_end`/`region_subseq`) ([#38](https://github.com/seanome/kmerseek/pull/38))
- Fixed amino acid ambiguity codes (B/J/Z) to resolve deterministically instead of by random draw, making indexing reproducible; U and O now map to C and K instead of being treated as unknown ([#39](https://github.com/seanome/kmerseek/pull/39))
- Built the k-mer frequency spectrum once instead of scanning it 9 times during indexing (up to 4.4x faster), and fixed gzip output silently truncating on write errors ([#40](https://github.com/seanome/kmerseek/pull/40))
- Fixed a panic when indexing with `--ksize 0`; k-mer sizes of 0 or above 100 are now rejected with a clear error instead of crashing or silently doing nothing ([#41](https://github.com/seanome/kmerseek/pull/41))
- Added `hp_lehninger_hpc`, a 3-letter hydrophobic/polar/cystine alphabet that gives cysteine its own symbol, for detecting C2H2 zinc finger proteins ([#42](https://github.com/seanome/kmerseek/pull/42))
