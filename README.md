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
and for the HP-family alphabets the **HP-encoded** window as well. The second catches
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
database (`hp_lehninger2`, k=17, every hit shown):

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

Use a `--ksize` of at least 12. Shorter k-mers under a reduced alphabet produce so many
overlapping sliding-window matches that the hit track becomes a wall of fragments.

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

The two- and three-class alphabets share an `hp_` prefix. They differ only on borderline
residues, so grouping them keeps a sweep easy to write and to grep. The multi-letter ones
use the names their papers use, digit included, so `sdm12` and `gbmr4` can be looked up.

Seeded negative controls put the seed after the class count, so `hp_random_control2_3`
is seed 3 of a 2-class control rather than a 23-class alphabet. `--random-seed 1..10`
gives independent replicates.

### Indexes built before this change

Only the current names parse. An index recording `protein`, `dayhoff`, `hp`, or an
`hp_<name>` without its class count fails to open:

```
Unknown alphabet in database: hp
```

The k-mers inside are still valid. Every alphabet hashes as it did before: `protein20`
and `dayhoff6` use the same sourmash hash functions as `protein` and `dayhoff`,
`hp_lehninger2` uses the same one as `hp`, and each `hp_<name>2` uses the same table as
`hp_<name>`. Only the recorded name changed, so rebuilding produces the same hits.

## HP Alphabet Variants

`--alphabet hp_lehninger2` collapses the 20 canonical amino acids down to hydrophobic (`h`)
/ polar (`p`) before k-mer extraction. The alphabets below (see
`src/rust/hp_alphabets.rs`) all agree on 15 of the 20 residues and differ only on
the five borderline ones: C, G, P, W and Y. `hp_lehninger2` is the one sourmash's own HP
encoding uses; the others follow the same `hp_<name>2` pattern and exist for the alphabet
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
and `DEHKNQRST` always polar. The schemes differ only in where they put C, G, P, W and Y.

`hp_lehninger_hpc3` has three classes rather than two. It keeps Lehninger's split for the
other 19 residues and gives cysteine its own symbol `c`.

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

Class count and k-size trade off against each other, so a k that suits a 2-class
alphabet is usually too long here. Searching CED-9 against the 25-sequence BCL-2 test
file at k=10 finds the same 21 targets under `hp_lehninger2` and `gbmr4`, but
`hp_lehninger2` reports 1673 matched regions against GBMR4's 294. `sdm12` finds nothing
at k=10, and 17 targets in 110 regions at k=5.

### Ambiguity codes

B (Asx), J (Xle) and Z (Glx) each stand for two residues. Every k-mer covering one is
indexed under both readings, so a query holding either residue matches. kmerseek does not
pick a representative: under `sdm12` and `hsdm17`, Asp and Asn fall in different classes,
so choosing one would assert a residue the source never had.

Whether that adds k-mers depends on the alphabet, because the two readings do not always
encode differently. Under `dayhoff6` Asp and Asn are both class `c`, and under the HP
tables both are polar, so the two readings of a window produce the same encoded k-mer and
the same hash. Nothing is added. Under `protein20`, `sdm12` and `hsdm17` they encode
differently, so each affected window contributes two k-mers instead of one.

`PLANTANDANIMALGENBMES` is 21 residues with a B at index 17, so at k=5 it has 17
windows and four of them span the B. Those four, with both readings and what each
encodes to:

| window | residues | readings | `dayhoff6` | `hp_lehninger2` |
|---|---|---|---|---|
| 13 | `LGENB` | `LGEND` `LGENN` | `ebccc` `ebccc` | `hhppp` `hhppp` |
| 14 | `GENBM` | `GENDM` `GENNM` | `bccce` `bccce` | `hppph` `hppph` |
| 15 | `ENBME` | `ENDME` `ENNME` | `cccec` `cccec` | `ppphp` `ppphp` |
| 16 | `NBMES` | `NDMES` `NNMES` | `ccecb` `ccecb` | `pphpp` `pphpp` |

Under `dayhoff6` and `hp_lehninger2` the two readings encode to the same string, so
the window still yields one k-mer. Under `protein20`, `sdm12` and `hsdm17` they encode
differently, so the window yields two k-mers instead of one.

| alphabet | distinct k-mers |
|---|---|
| `protein20` | 21 |
| `sdm12` | 21 |
| `hsdm17` | 21 |
| `dayhoff6` | 17 |
| `hp_lehninger2` | 14 |

Where each count comes from:

- `protein20` has 17 k-mers before expansion, one per window, all distinct. Expanding
  B into D and N turns each of the four B-windows into two k-mers, adding 4. That gives
  21 k-mers.
- `sdm12` and `hsdm17` behave the same way, 17 k-mers before expansion and 21 after,
  because they also keep Asp and Asn in separate classes.
- `dayhoff6` has 17 k-mers before expansion. Expanding B adds none, because `LGEND` and
  `LGENN` both encode to `ebccc`. It stays at 17 k-mers.
- `hp_lehninger2` has 14 k-mers before expansion, not 17. With only two symbols, three
  pairs of ordinary windows already encode identically: `PLANT` and `ALGEN` are both
  `hhhpp`, `ANTAN` and `ANDAN` are both `hpphp`, `NTAND` and `NBMES` are both `pphpp`.
  Expanding B adds none, for the same reason as `dayhoff6`. It stays at 14 k-mers.

So the only alphabets where de-ambiguating B changes the k-mer count are the ones that
encode D and N differently. The gap between 17 and 14 is HP collapsing distinct
residues, which would be there with or without the B.

An ambiguity code doubles the readings of every window it falls in, so a window holding
*n* codes expands to 2^*n*. kmerseek indexes a window holding at most 4 codes, which is
16 readings. A window with more is dropped: indexing only part of its readings would make
matching depend on which subset was kept, which is worse than not indexing that one
window. SwissProt holds about 900 non-canonical residues in 207.6 M, so a window with
five codes should not arise.

U (Sec) and O (Pyl) are handled differently. They are specific residues rather than
ambiguities, so each takes its closest canonical analogue, C and K, under a reduced
alphabet. Under `protein20` they are kept as themselves.

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
