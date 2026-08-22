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
kmerseek index -i proteome.fasta --ksize 10 --alphabet hp_lehninger2 --remove-low-complexity
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
    --alphabet hp_lehninger2 --ksize 17

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

## Alphabets

Pick one with `--alphabet` (`-a`). Every name ends in the number of classes it collapses
the 20 amino acids into, so a filename or a results CSV says how much reduction happened:

| Alphabet | Classes | Source |
|---|:---:|---|
| `protein20` | 20 | the full alphabet, no reduction |
| `dayhoff6` | 6 | Dayhoff |
| `hp_lehninger2` | 2 | Lehninger, the HP default |
| `hp_thomas_dill2` | 2 | Thomas & Dill 1996 |
| `hp_kyte_doolittle2` | 2 | Kyte & Doolittle 1982 |
| `hp_thomas_dill_no_c2` | 2 | Thomas-Dill, C reassigned to polar |
| `hp_lehninger_c_nonpolar2` | 2 | Lehninger, C reassigned to hydrophobic |
| `hp_lehninger_hpc3` | 3 | Lehninger, C given its own class |
| `hp_pbotc_1st_ed2` | 2 | Physical Biology of the Cell, 1st ed |
| `hp_random_control2` | 2 | negative control, randomized h/p split |
| `gbmr4` | 4 | Solis & Rackovsky 2000 |
| `wwmj5` | 5 | Wang & Wang 1999 |
| `gbmr7` | 7 | Solis & Rackovsky 2000 |
| `sdm12` | 12 | Prlic et al. 2000 |
| `mmseqs12` | 12 | Steinegger & Soding 2018 |
| `wass14` | 14 | Ieremie et al. 2024 |
| `hsdm17` | 17 | Prlic et al. 2000 |
| `uniprot18` | 18 | Ieremie et al. 2024 |

The two-and-three-class alphabets keep an `hp_` prefix because they are one family that
differs only on borderline residues, and grouping them makes a sweep easy to write. The
multi-letter ones use the names their papers use, digit included, so `sdm12` and `gbmr4`
are directly citeable.

Seeded negative controls put the seed after the class count:
`hp_random_control2_3` is seed 3 of a 2-class control, not a 23-class alphabet. Use
`--random-seed 1..10` to get independent replicates.

### Older spellings

These all still work on the command line and in existing indexes, and are normalized to
the current name on read:

| Older | Now |
|---|---|
| `protein`, `raw` | `protein20` |
| `dayhoff` | `dayhoff6` |
| `hp_<name>` without a class count | `hp_<name>2` |
| a `reduced_` prefix on anything | the same name without it |
| `shuffled_control` | `random_control` |
| `--encoding` | `--alphabet` |
| `--shuffled-seed` | `--random-seed` |

Nothing writes the older names any more.

### Indexes built before the rename

Most keep working. `dayhoff`, `protein` and the `hp_<name>` alphabets resolve to the same
hash function and the same table as before, so their k-mer hashes are unchanged and
existing databases are searchable as-is.

The one exception is `hp`. It used to select sourmash's built-in HP encoder, which applies
the same Lehninger partition but hashes the encoded sequence as lowercase `h`/`p`, where
kmerseek's own tables are uppercased to `H`/`P` first. The two share no k-mer hashes at
all. `hp` is now a spelling of `hp_lehninger2`, so an index built under the old `hp` would
sketch queries on the uppercase path and match nothing. kmerseek refuses to open such an
index and tells you to rebuild it:

```
This index was built with the old `hp` encoding, whose k-mer hashes are not compatible
with `hp_lehninger2` ... Rebuild the index with --alphabet hp_lehninger2
```

Rebuilding reproduces the same hits: the amino-acid partition never changed, only the
bytes that get hashed.

## HP Alphabet Variants

`--alphabet hp_lehninger2` collapses the 20 canonical amino acids down to hydrophobic (`h`)
/ polar (`p`) before k-mer extraction. The alphabets below (see
`src/rust/hp_alphabets.rs`) all agree on 15 of the 20 residues and differ only on
the five borderline ones: C, G, P, W and Y. Lehninger is the current
default, spelled `hp_lehninger2`; the others are selectable via the same
`hp_<name>2` pattern (e.g. `hp_thomas_dill2`) for the alphabet
robustness sweep.

Cysteine is the residue the schemes disagree about most, because its thiol side chain is
nonpolar but its disulfide bonding is a chemistry of its own. `hp_lehninger_c_nonpolar2`
folds it into `h`; `hp_lehninger_hpc3` gives it a third class instead.

Alphabet / Scheme | Hydrophobic (`h`) | Polar (`p`)
-- | -- | --
Lehninger (default) | `AFILMV` GPWY | `DEHKNQRST` C
Thomas-Dill / PBotC 2nd | `AFILMV` CWY | `DEHKNQRST` GP
Kyte-Doolittle | `AFILMV` C | `DEHKNQRST` GPWY
TD−C | `AFILMV` WY | `DEHKNQRST` CGP
Leh+C | `AFILMV` CGPWY | `DEHKNQRST`
PBotC 1st | `AFILMV` CPWY | `DEHKNQRST` G

The residues in backticks are fixed across every scheme: `AFILMV` is always hydrophobic
and `DEHKNQRST` always polar. Only the five borderline residues C, G, P, W and Y move,
which is the whole of the disagreement between these alphabets.

`hp_lehninger_hpc3` is the odd one out, with three classes rather than two: it keeps
Lehninger's split for the other 19 residues and gives cysteine its own symbol `c`.

| Alphabet / Scheme | Hydrophobic (`h`) | Polar (`p`) | Cystine (`c`) |
| -- | -- | -- | -- |
| Leh HPC | `AFILMV` GPWY | `DEHKNQRST` | C |

## Multi-Letter Reduced Alphabets

The HP alphabets above answer one question per residue. The alphabets in this section
(see `src/rust/reduced_alphabets.rs`) keep 4 to 18 classes, so they discard less
chemistry per position while still collapsing the substitutions that proteins tolerate
most often.

Peterson et al. (2009) benchmarked over 150 published clustering schemes against DALI
fold assignments and found that reduced alphabets beat the full 20-letter alphabet,
with the best results at 9-12 classes. Ieremie et al. (2024) reused their top three and
added five more when testing how alphabet reduction affects protein language models.
Where the two papers overlap, their partitions are identical.

| Moltype | Classes | Clusters | Source |
|---------|:---:|---|---|
| `gbmr4` | 4 | `ADKERNTSQ` `YFLIVMCWH` `G` `P` | Solis & Rackovsky 2000 |
| `wwmj5` | 5 | `CMFILVWY` `ATH` `GP` `DE` `SNQRK` | Wang & Wang 1999 |
| `gbmr7` | 7 | `DN` `AEFIKLMQRVWY` `CH` `T` `S` `G` `P` | Solis & Rackovsky 2000 |
| `sdm12` | 12 | `A` `D` `KER` `N` `TSQ` `YF` `LIVM` `C` `W` `H` `G` `P` | Prlic et al. 2000 |
| `mmseqs12` | 12 | `AST` `LM` `IV` `KR` `EQ` `ND` `FY` `C` `G` `H` `P` `W` | Steinegger & Soding 2018 |
| `wass14` | 14 | `WM` `DI` `P` `C` `AV` `K` `T` `RE` `G` `L` `Y` `SH` `F` `NQ` | Ieremie et al. 2024 |
| `hsdm17` | 17 | `A` `D` `KE` `R` `N` `T` `S` `Q` `Y` `F` `LIV` `M` `C` `W` `H` `G` `P` | Prlic et al. 2000 |
| `uniprot18` | 18 | `A` `R` `N` `D` `C` `Q` `EP` `G` `HL` `I` `K` `M` `F` `S` `T` `W` `Y` `V` | Ieremie et al. 2024 |

GBMR4, SDM12 and HSDM17 were the top performers in Peterson et al. on recall at 0.01
errors per query, AUC and mean pooled precision respectively. Each is a refinement of
the previous one: going from GBMR4 to SDM12 to HSDM17 only splits classes, never moves
a residue across an existing boundary (`test_hsdm17_refines_sdm12_refines_gbmr4`).

Each class is written as its first residue in lowercase, so an SDM12-encoded sequence
shows `LIVM` as `l` and can be read against the source residues directly.

Class count and k-size trade off against each other, so a k that works for `hp` is
usually too long here. Searching CED-9 against the 25-sequence BCL-2 test file at k=10
finds the same 21 targets under `hp` and `gbmr4`, but `hp` reports 1673 matched
regions against GBMR4's 294; `sdm12` finds nothing at k=10 and 17 targets in 110
regions at k=5.

### Ambiguity codes

Under `dayhoff6` and the named HP tables, B (Asx), J (Xle) and Z (Glx) are replaced
by a fixed representative, because both residues each code stands for land in the same
class either way. Most of these alphabets split at least one of those pairs -- SDM12 and
HSDM17 give Asp and Asn separate classes -- so under them B/J/Z, along with U (Sec) and
O (Pyl), are kept and hashed as themselves, the way X already is. Only `gbmr4`
and `gbmr7` collapse all three pairs and so still substitute.

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
