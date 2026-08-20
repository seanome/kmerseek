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

## HP Alphabet Variants

`--encoding hp` collapses the 20 canonical amino acids down to hydrophobic (`h`)
/ polar (`p`) before k-mer extraction. The alphabets below (see
`src/rust/hp_alphabets.rs`) all agree on 15 of the 20 residues and differ only on
the five borderline ones -- **C, G, P, W, Y** (bolded). Lehninger is the current
default (`hp` moltype); the others are selectable via `hp_<name>` moltypes
(e.g. `hp_thomas_dill`) for the alphabet robustness sweep.

`hp_lehninger_hpc` is a 3-letter variant: it keeps Lehninger's H/P split for
every residue except cysteine, which gets its own third symbol `c` (cystine)
instead of being folded into `h` the way `hp_lehninger_c_nonpolar` does --
disulfide-bond formation is a distinct chemistry from ordinary hydrophobic
packing.

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
