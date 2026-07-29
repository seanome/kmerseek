# Changelog

## 0.4.0

- Added a simple visualization of target hits against a query sequence, showing both the original and encoded alphabets ([#28](https://github.com/seanome/kmerseek/pull/28))
- Added removal of hits with uncorrected p-value > 0.05, since those are removed anyway after BH correction ([#27](https://github.com/seanome/kmerseek/pull/27))
- Added filtering by at least 2 shared hashes per hit to remove spurious hits ([#27](https://github.com/seanome/kmerseek/pull/27))
- Added a table of HP encodings to the README ([#33](https://github.com/seanome/kmerseek/pull/33))
- Removed the `--scaled` option since we always use `--scaled 1` anyway ([#30](https://github.com/seanome/kmerseek/pull/30))
- Sorted top hits in the visualization by Benjamini-Hochberg corrected p-value (q-value) instead of containment ([#32](https://github.com/seanome/kmerseek/pull/32))
