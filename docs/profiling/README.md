# Where a search spends its time (Swiss-Prot, October 2026)

![Seconds per search step, protein20 k=10 and hp_lehninger2 k=19 on Swiss-Prot](../images/search_step_profile_swissprot.png)

100 human proteins searched against all of Swiss-Prot (575,748 proteins, 209 million
residues), on one Sherlock node with 8 threads, at commit 07b7244. Every step's seconds
come from the timers kmerseek turns on when `KMERSEEK_PROFILE` is set. Steps that run on
8 threads at once (orange) are summed over the threads, so they add up to more than the
wall time. Bars are the mean of two runs with the timers on.

| Index | Wall time | Peak memory | Candidate targets per query | Results kept |
|---|---|---|---|---|
| protein20, k=10 | 105 s | 15 GB | 27 | 1,712 |
| hp_lehninger2, k=19 | 181 s | 55 GB | 66,658 (11.6% of Swiss-Prot) | 4,065,945 |

## protein20, k=10: loading the index is 72% of the run

Loading the index takes 75 s of the 105 s wall time. The index holds 117 million different
k-mers, and `load_search_cache` reads every one of them into an in-memory `HashMap` before
the first query starts. The 100 queries themselves take 3.3 thread-seconds, about 0.4 s of
wall time on 8 threads. The search would be about 3 times faster with an index that is
read from disk as needed (sorted k-mers in a memory-mapped file) instead of rebuilt in memory
on every run.

## hp_lehninger2, k=19: too many candidates

With two letters and k=19 there are only 2^19 = 524,288 possible k-mers. The Swiss-Prot
index holds 572,309 different hashes; the extra ones presumably contain letters outside the
two classes, which was not checked. Each query k-mer therefore points at many
proteins: one query has on average 66,658 candidate targets, and 61% of those pass the
default filters and become results (6.5 million CSV rows). The time follows the number of
candidates:

- finding matched regions: 606 thread-seconds, 58% of the search work;
- the shared k-mers between query and target: 292 thread-seconds, 28%;
- reading targets from RocksDB: 130 thread-seconds, 12%. 446,070 different targets were
  read, 77% of Swiss-Prot. All of them stay in the in-memory cache, which is why the run
  peaks at 55 GB.

Finding the candidates in the inverted index takes only 0.6 s. What would help is fewer
candidates reaching the per-target steps: a cap on k-mers that occur in too many proteins,
or requiring two shared k-mers on the same diagonal (target position minus query position)
before a target is read from disk.

## What the timers do not cover

About 29 s of each run's wall time is outside every timer, in both settings. That is the
wall time minus loading, writing the CSV, and the parallel steps divided by 8. Freeing the
index and cache at exit is a likely place, but this was not measured.

The timers themselves cost nothing measurable: runs with them off took 104 s and 191 s,
runs with them on 105-107 s and 175-177 s. The first run of each setting (140 s and 466 s)
read the index from a cold disk cache and is left out of the means.

## Files

- `timing.tsv`: wall time, CPU time, peak memory and CSV rows for every run.
- `steps.tsv`: seconds and number of calls per step, for the runs with timers on.
- `commit.txt`, `inputs.md5`: the kmerseek commit and the md5 of Swiss-Prot and the queries.
- Rerun on Sherlock: `sbatch scripts/profile_search_steps.sbatch`; redraw:
  `python scripts/plot_search_profile.py`.

# Filtering candidates before reading targets (hp_lehninger2, k=19)

![Thread-seconds per step before and after counting shared k-mers first, and wall time per k-mer cap](../images/candidate_filter_swissprot_hp_lehninger2_k19.png)

Commit fa352ce counts each candidate's shared k-mers from the posting lists. It drops
candidates below `--min-shared-kmers` (2 by default) or the containment threshold before
their signature is read from disk. A posting list names each target once, so the count is
exactly the number `compare` checks later. The search results do not change: the CSVs are
byte-identical to the runs above (`candidate_filter/identical.md5`).

- Targets compared per query fell from 66,658 to 40,659. The 39% dropped shared only one k-mer.
- Time spent on shared k-mers fell from 292 to 243 thread-seconds, and reading targets from
  130 to 102. Target reads fell only from 446,070 to 435,401, because most dropped targets
  are still read for another query.
- Finding matched regions is unchanged at 604 thread-seconds: every pair that reaches it
  was going to be a result.
- Wall time went from 181 s to 167 s (8% faster); CPU time from 1,012 s to 952 s.

`--max-kmer-targets N` (off by default) keeps a pair only if it shares at least one k-mer
found in no more than N targets. On this index it does not pay:

| Cap | Wall time | Pairs kept | Top 1% of pairs by region score kept |
|---|---|---|---|
| none | 167 s | 4,065,945 (100%) | 100% |
| 10,000 | 171 s | 99.7% | 100% |
| 1,000 | 160 s | 89.0% | 98.8% |
| 100 | 2 s | 0 | 0 |

A cap of 100 leaves nothing because every k-mer of these 100 queries occurs in more than 100
Swiss-Prot proteins. With 2^19 possible k-mers and 209 million residues, no hp_lehninger2
k=19 k-mer is rare. A cap of 1,000 costs 11% of the pairs for 4% less time.

The work is per pair, and the 4 million pairs pass the default filter of two shared k-mers.
In an alphabet this small, two shared k-mers are expected by chance alone. Two changes could
cut the pairs:

- a default minimum that grows with the number of k-mers two random proteins would share;
- storing each k-mer's position in the posting lists, so pairs whose shared k-mers do not lie
  on one diagonal can be dropped before the target is read. This changes the index format.

Both change which results are reported, so neither is in this PR. Rerun with
`sbatch scripts/profile_candidate_filter.sbatch` and redraw with
`python scripts/plot_candidate_filter.py`.

## Posting files read in place

![Search wall time, peak memory and index load time, before and after](../images/posting_files_swissprot.png)

Before this change, search read every k-mer of the index into three HashMaps before the
first query ran: the posting lists, each k-mer's frequency, and each k-mer's IDF. For
protein20 k=10 that is 116,913,137 k-mers; for hp_lehninger2 k=19, 572,309. The index now
keeps its posting lists in four flat files under `<index>/postings/`, and search maps them
into memory instead of reading them (`src/rust/postings.rs`). A lookup goes to the k-mer's
shard and binary-searches it. Opening checks every shard and offset boundary, which reads the
whole offsets file (0.94 GB for protein20) but none of the k-mers or targets.

Sherlock job 50498949 ran the old binary on the old indexes and the new binary (commit
3cddf3b) on rebuilt indexes, alternating, three runs each, on one node. The same 100 human
proteins were searched against Swiss-Prot. Run 1 of each read the index from a cold disk
cache.

| | protein20 k=10, before | after | hp_lehninger2 k=19, before | after |
|---|---|---|---|---|
| Wall time, warm runs 2 and 3 | 103 s, 109 s | 1.0 s, 1.0 s | 165 s, 184 s | 185 s, 175 s |
| Wall time, run 1 (cold) | 240 s | 13 s | 456 s | 421 s |
| CPU time, warm runs (mean) | 107 s | 2.4 s | 957 s | 1,048 s |
| Peak memory | 15.0 GB | 2.0 GB | 55.0 GB | 53.8 GB |
| Load the index | 74-206 s | 0.4-1.0 s | 0.9-49 s | 0.05-0.7 s |
| Index on disk | 6.29 GB | 7.04 GB | 4.88 GB | 4.88 GB |

All twelve result CSVs match their "before" run byte for byte (`posting_files/identical.md5`):
7,475 rows for protein20, 6,542,668 for hp_lehninger2.

- protein20 k=10 was slow because of loading the index, so it gains the most: about 100
  times faster warm, and 7 times less memory.
- hp_lehninger2 k=19 loaded its index in under a second already. Its time goes to the
  4,065,945 query-target pairs (PR 158), so its wall time does not change beyond the
  run-to-run noise. It uses 10% more CPU: the "shared k-mers" step went from 243 to 334
  thread-seconds (mean of the warm runs). That step looks up the database frequency of every
  shared k-mer twice per pair (`calculate_mean_matched_kmer_freq` and
  `calculate_sum_matched_kmer_freq`), and each lookup is now a binary search within a shard
  of about 560 k-mers instead of one hash lookup. Its 54 GB peak is the cache of target
  sketches read from disk, which this change does not touch.
- protein20's peak memory after the change is 2.0 GB here and was 1.2 GB in job 48806130,
  before the boundary check was added in review. The difference is close to the size of the
  offsets file the check reads (0.94 GB); that cause was not tested separately.
- The protein20 index grew by 0.75 GB. Its posting files are 2.60 GB, stored uncompressed;
  the old `ii_shard_{s}` values were RocksDB values, which RocksDB compresses by default. The
  hp_lehninger2 index is the same size before and after; its posting files are 0.68 GB.
- Each index rebuilt in 2 min 14 s with a peak under 1.4 GB.
- The "before" warm runs vary between jobs (protein20: 105 s in job 48672829, 165-188 s in
  job 48806130, 103-109 s here), with the same binary and index. Both versions run in the
  same job, so the comparison within a job holds.

Rerun with `sbatch scripts/profile_posting_files.sbatch` and redraw with
`python scripts/plot_posting_files.py`.
