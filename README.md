# Kmerseek

## Compiling on Mac

You may need to add these magical `export` commands to make your Python install work:

```bash
export MACOSX_DEPLOYMENT_TARGET=10.15 \
	&& export PYTHON_CONFIGURE_OPTS="--enable-framework" \
	&& export PYTHON_SYS_EXECUTABLE="$(which python)" \
	&& export PYO3_PYTHON="$(which python)" \
	&& export PYTHONPATH="/Users/olga/anaconda3/envs/kmerseek-dev/lib/python3.13/site-packages:$PYTHONPATH" \
	&& export DYLD_FALLBACK_LIBRARY_PATH="/Users/olga/anaconda3/envs/kmerseek-dev/lib:$DYLD_FALLBACK_LIBRARY_PATH" \
	&& export RUSTFLAGS="-C link-arg=-undefined -C link-arg=dynamic_lookup"
```

It may look like this:

```bash
export MACOSX_DEPLOYMENT_TARGET=10.15 \
        && export PYTHON_CONFIGURE_OPTS="--enable-framework" \
        && export PYTHON_SYS_EXECUTABLE="/Users/olga/anaconda3/envs/kmerseek-dev/bin/python" \
        && export PYO3_PYTHON="/Users/olga/anaconda3/envs/kmerseek-dev/bin/python" \
        && export PYTHONPATH="/Users/olga/anaconda3/envs/kmerseek-dev/lib/python3.13/site-packages:$PYTHONPATH" \
        && export DYLD_FALLBACK_LIBRARY_PATH="/Users/olga/anaconda3/envs/kmerseek-dev/lib:$DYLD_FALLBACK_LIBRARY_PATH" \
        && export RUSTFLAGS="-C link-arg=-undefined -C link-arg=dynamic_lookup"
```

## Testing

Run real-world examples like:

```bash
cargo run --example test_bcl2_processing
```

## Removing low-complexity k-mers

Low-complexity k-mers -- homopolymer runs like a poly-glutamate tract (`EEEEE`) or,
under a reduced alphabet, an all-hydrophobic window (`hhhhh`) -- are abundant and
carry little discriminative signal. Pass `--remove-low-complexity` at index time to
drop them:

```bash
kmerseek index -i proteome.fasta --ksize 10 --encoding hp --remove-low-complexity
```

Two independent checks run per k-mer: the **raw amino-acid** window (any encoding),
and for `hp`-family encodings the **HP-encoded** window as well. The second catches
windows that aren't raw homopolymers but still collapse to one symbol -- `LIVMA` is
five different residues that all encode to `h`.

Indexing reports what was removed, so you can tell whether the flag mattered:

```
Removed 75 of 9063 k-mer windows as low-complexity (0.83%)
```

The setting is **stored in the index**, and `kmerseek search` reads it back and
builds query sketches the same way. You don't repeat the flag when searching, and
search says plainly what the index holds and what it is doing:

```
Index: low-complexity k-mers were REMOVED when it was built
This search: low-complexity k-mers are REMOVED from query sketches (matching the index)
```

Pass `--remove-low-complexity` (or `--remove-low-complexity false`) to `search`
only to override that deliberately. Disagreeing with the index is allowed but
warned about:

```
This search: low-complexity k-mers are REMOVED from query sketches (--remove-low-complexity)
WARNING: this disagrees with the index. Containment is intersection / query_size,
so k-mers present on only one side still count toward the denominator and skew scores.
```

Results carry the setting too, in a `remove_low_complexity` column next to
`ksize`/`scaled`/`moltype`, so a CSV is self-describing without the command line
that produced it. It is a column rather than a `#` comment line because a comment
would break `pl.scan_csv` and every other plain CSV reader.

That symmetry matters. Containment is `intersection / query_size`, so if the index
dropped these k-mers but queries kept them, they'd match nothing while still
inflating the denominator -- deflating scores for exactly the queries containing
low-complexity regions.

Auto-generated filenames gain a `.nolowcomplexity` segment, so builds with and
without removal coexist instead of overwriting each other:

```
proteome.fasta.hp.k10.scaled1.kmerseek.rocksdb                  # default
proteome.fasta.hp.k10.scaled1.nolowcomplexity.kmerseek.rocksdb  # --remove-low-complexity
```

Removal is **off by default**; existing indexes and workflows are unaffected.
Note that only *exact* homopolymers are dropped -- a near-homopolymer such as
`hhhhhhhhhp` is kept.

## Visualizing hits

`scripts/visualize_hits.py` renders a per-gene PNG+SVG pair showing every hit
mapped onto the query protein: a full-length bar, each matched target's *actual*
matched regions positioned to scale (stacked into lanes when hits overlap, numbered
inside each box -- never floating text that can collide with a neighbor), and the
query / encoded-alphabet / target alignment printed beneath each hit, with every
region of a multi-region hit shown individually (not just one representative).
Each hit's containment, Jaccard, fold-enrichment, Poisson p-value, and a
Benjamini-Hochberg FDR-corrected q-value (corrected across every target actually
tested for that query, not just the ones kept after `--min-containment` filtering)
are printed alongside it. Categorical colors always come from a built-in matplotlib qualitative colormap,
sized to how many distinct targets there are (Set2 for up to 8, tab10 up to 10,
Set3 up to 12, tab20 beyond that) -- the alignment block's title is colored to
match its box. Example, `ced9.fasta` searched against a 25-protein BCL2-family
database (`hp` encoding, k=17, scaled=1, every hit shown):

![Example hit visualization for CED9_CAEEL](docs/images/ced9_hits_example.png)

([SVG version](docs/images/ced9_hits_example.svg) -- sequence text stays selectable/copyable)

```bash
kmerseek search -q ced9.fasta -t bcl2_family.rocksdb -o results.csv \
    --encoding hp --ksize 17 --scaled 1

python scripts/visualize_hits.py \
    --csv results.csv \
    --query-fasta ced9.fasta \
    --output-dir hits_png/
```

Use **`--scaled 1`** -- higher scaled values subsample k-mers and can silently drop
real hits (e.g. CED9 vs. BCL2_HUMAN itself disappears at `--scaled 2`, even though
it's a genuine match at `--scaled 1`). Use `--ksize` of at least 12 too -- shorter
k-mers under reduced alphabets like `hp` produce so many overlapping sliding-window
matches that the hit track becomes an unreadable wall of fragments.

`--query-fasta` supplies the full-length protein for the top bar and must be the
same FASTA used as the search query. Omit `--query-name` to render one PNG+SVG pair
per query found in the CSV. `--min-containment` and `--max-hits` are off by default
(every hit, for every target, is drawn); `--max-hits N` caps the figure to the top N
*distinct targets* by BH-corrected q-value, most significant first (all of a kept
target's hit spans are still shown, so one heavily-fragmented target can't crowd out
the others) -- use it to tame proteome-scale searches where a gene can have dozens of
distinct hits. See
`python scripts/visualize_hits.py --help` for all options.

## Encoding Names

Every reduced alphabet's moltype ends in the number of classes it collapses the 20
amino acids into, so the name states how much chemistry it discards:

| Moltype | Classes | Previously |
|---------|:---:|---|
| `reduced_dayhoff6` | 6 | `dayhoff` |
| `reduced_hp_<name>2` / `3` | 2 or 3 | `hp_<name>` |
| `reduced_gbmr4` ... `reduced_uniprot18` | 4-18 | new |

The old spellings are still accepted on the command line and in existing indexes, and
are normalized to the current name on read, so databases built before the rename keep
working. Nothing writes the old names any more.

`protein` (the full 20-letter alphabet) and `hp` are unchanged. `hp` is sourmash's
built-in HP encoding: it uses the same Lehninger partition as
`reduced_hp_lehninger2`, but sourmash hashes it as lowercase `h`/`p` while the custom
tables are uppercased to `H`/`P` before hashing, so the two share no k-mer hashes and
are not interchangeable. Pick one and stay with it for a given index.

## HP Alphabet Variants

`--encoding hp` collapses the 20 canonical amino acids down to hydrophobic (`h`)
/ polar (`p`) before k-mer extraction. The alphabets below (see
`src/rust/hp_alphabets.rs`) all agree on 15 of the 20 residues and differ only on
the five borderline ones -- **C, G, P, W, Y** (bolded). Lehninger is the current
default (`hp` moltype); the others are selectable via `reduced_hp_<name>2` moltypes
(e.g. `reduced_hp_thomas_dill2`) for the alphabet robustness sweep.

The trailing digit is the class count, matching the multi-letter alphabets in the
next section, so every moltype states its size. These were previously named
`hp_<name>` without the count. The old spellings are still accepted on the command
line and in existing indexes, and are normalized to the current name on read, so a
database built before the rename keeps working; nothing writes them any more.

`reduced_hp_lehninger_hpc3` is a 3-letter variant: it keeps Lehninger's H/P split for
every residue except cysteine, which gets its own third symbol `c` (cystine)
instead of being folded into `h` the way `reduced_hp_lehninger_c_nonpolar2` does --
disulfide-bond formation is a distinct chemistry from ordinary hydrophobic
packing.

| Moltype | Classes | Previously |
|---------|:---:|---|
| `reduced_hp_lehninger2` | 2 | `hp_lehninger` |
| `reduced_hp_thomas_dill2` | 2 | `hp_thomas_dill` |
| `reduced_hp_kyte_doolittle2` | 2 | `hp_kyte_doolittle` |
| `reduced_hp_thomas_dill_no_c2` | 2 | `hp_thomas_dill_no_c` |
| `reduced_hp_lehninger_c_nonpolar2` | 2 | `hp_lehninger_c_nonpolar` |
| `reduced_hp_lehninger_hpc3` | 3 | `hp_lehninger_hpc` |
| `reduced_hp_pbotc_1st_ed2` | 2 | `hp_pbotc_1st_ed` |
| `reduced_hp_shuffled_control2` | 2 | `hp_shuffled_control` |

Seeded controls put the seed after the class count -- `reduced_hp_shuffled_control2_3`
is seed 3 of a 2-class control, not a 23-class alphabet.

| AA | Lehninger (current) | Thomas-Dill/PBotC 2nd | Kyte-Doolittle | TD−C | Leh+C | Leh HPC (3-letter) | PBotC 1st |
|----|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| A | h | h | h | h | h | h | h |
| **C** | p | h | h | p | h | c | h |
| D | p | p | p | p | p | p | p |
| E | p | p | p | p | p | p | p |
| F | h | h | h | h | h | h | h |
| **G** | h | p | p | p | h | h | p |
| H | p | p | p | p | p | p | p |
| I | h | h | h | h | h | h | h |
| K | p | p | p | p | p | p | p |
| L | h | h | h | h | h | h | h |
| M | h | h | h | h | h | h | h |
| N | p | p | p | p | p | p | p |
| **P** | h | p | p | p | h | h | h |
| Q | p | p | p | p | p | p | p |
| R | p | p | p | p | p | p | p |
| S | p | p | p | p | p | p | p |
| T | p | p | p | p | p | p | p |
| V | h | h | h | h | h | h | h |
| **W** | h | h | p | h | h | h | h |
| **Y** | h | h | p | h | h | h | h |

## Multi-Letter Reduced Alphabets

The HP alphabets above answer one question per residue. The alphabets in this section
(see `src/rust/reduced_alphabets.rs`) keep 4 to 18 classes, so they discard less
chemistry per position while still collapsing the substitutions that proteins tolerate
most often. They follow the same naming rule: the number in the moltype is the class
count.

Peterson et al. (2009) benchmarked over 150 published clustering schemes against DALI
fold assignments and found that reduced alphabets beat the full 20-letter alphabet,
with the best results at 9-12 classes. Ieremie et al. (2024) reused their top three and
added five more when testing how alphabet reduction affects protein language models.
Where the two papers overlap, their partitions are identical.

| Moltype | Classes | Clusters | Source |
|---------|:---:|---|---|
| `reduced_gbmr4` | 4 | `ADKERNTSQ` `YFLIVMCWH` `G` `P` | Solis & Rackovsky 2000 |
| `reduced_wwmj5` | 5 | `CMFILVWY` `ATH` `GP` `DE` `SNQRK` | Wang & Wang 1999 |
| `reduced_gbmr7` | 7 | `DN` `AEFIKLMQRVWY` `CH` `T` `S` `G` `P` | Solis & Rackovsky 2000 |
| `reduced_sdm12` | 12 | `A` `D` `KER` `N` `TSQ` `YF` `LIVM` `C` `W` `H` `G` `P` | Prlic et al. 2000 |
| `reduced_mmseqs12` | 12 | `AST` `LM` `IV` `KR` `EQ` `ND` `FY` `C` `G` `H` `P` `W` | Steinegger & Soding 2018 |
| `reduced_wass14` | 14 | `WM` `DI` `P` `C` `AV` `K` `T` `RE` `G` `L` `Y` `SH` `F` `NQ` | Ieremie et al. 2024 |
| `reduced_hsdm17` | 17 | `A` `D` `KE` `R` `N` `T` `S` `Q` `Y` `F` `LIV` `M` `C` `W` `H` `G` `P` | Prlic et al. 2000 |
| `reduced_uniprot18` | 18 | `A` `R` `N` `D` `C` `Q` `EP` `G` `HL` `I` `K` `M` `F` `S` `T` `W` `Y` `V` | Ieremie et al. 2024 |

GBMR4, SDM12 and HSDM17 were the top performers in Peterson et al. on recall at 0.01
errors per query, AUC and mean pooled precision respectively. Each is a refinement of
the previous one: going from GBMR4 to SDM12 to HSDM17 only splits classes, never moves
a residue across an existing boundary (`test_hsdm17_refines_sdm12_refines_gbmr4`).

Each class is written as its first residue in lowercase, so an SDM12-encoded sequence
shows `LIVM` as `l` and can be read against the source residues directly.

Class count and k-size trade off against each other, so a k that works for `hp` is
usually too long here. Searching CED-9 against the 25-sequence BCL-2 test file at k=10
finds the same 21 targets under `hp` and `reduced_gbmr4`, but `hp` reports 1673 matched
regions against GBMR4's 294; `reduced_sdm12` finds nothing at k=10 and 17 targets in 110
regions at k=5.

### Ambiguity codes

Under `reduced_dayhoff6`, `hp` and the named HP tables, B (Asx), J (Xle) and Z (Glx) are replaced
by a fixed representative, because both residues each code stands for land in the same
class either way. Most of these alphabets split at least one of those pairs -- SDM12 and
HSDM17 give Asp and Asn separate classes -- so under them B/J/Z, along with U (Sec) and
O (Pyl), are kept and hashed as themselves, the way X already is. Only `reduced_gbmr4`
and `reduced_gbmr7` collapse all three pairs and so still substitute.

## Using the Builder Pattern

The `ProteomeIndex` now supports a fluent Builder pattern:

```rust
use kmerseek::index::ProteomeIndex;

// Using the builder pattern
let index = ProteomeIndex::builder()
    .path("/path/to/database.db")
    .ksize(5)
    .scaled(1)
    .moltype("protein")
    .build()?;

// With auto filename generation
let index = ProteomeIndex::builder()
    .path("/path/to/base")
    .ksize(5)
    .scaled(1)
    .moltype("protein")
    .build_with_auto_filename()?;

// With raw sequence storage
let index = ProteomeIndex::builder()
    .path("/path/to/database.db")
    .ksize(5)
    .scaled(1)
    .moltype("protein")
    .store_raw_sequences(true)
    .build()?;

// Dropping low-complexity (homopolymer) k-mers
let index = ProteomeIndex::builder()
    .path("/path/to/database.db")
    .ksize(5)
    .scaled(1)
    .moltype("hp")
    .remove_low_complexity(true)
    .build()?;
```

You can also use convenience methods:

```rust
// Create a new index
let index = ProteomeIndex::new_simple(
    "/path/to/database.db",
    5,        // k-mer size
    1,        // scaled
    "protein", // molecular type
    false,    // don't store raw sequences
)?;

// With auto filename generation
let index = ProteomeIndex::new_with_auto_filename_simple(
    "/path/to/data.fasta",
    5,        // k-mer size
    1,        // scaled
    "protein", // molecular type
    false,    // don't store raw sequences
)?;
```

Run the builder pattern demo:

```bash
cargo run --example builder_pattern_demo
```
