# What Olga flags in review

Olga's own review comments on seanome/kmerseek PRs #36 to #75, gathered on 2026-10-01 from
`gh api repos/seanome/kmerseek/pulls/comments` (her comments only, coding-agent replies dropped).
Each rule below is something she had to ask for by hand at least once. A review that finds
none of these missing has checked for each one.

To refresh this list later:

```bash
GH_PAGER=cat gh api --paginate "repos/seanome/kmerseek/pulls/comments?per_page=100" \
  --jq '.[] | select(.user.login=="olgabot") | select(.body | test("Written by Claude|Claude Code|🤖") | not)
        | "PR\(.pull_request_url|split("/")|last) \(.path):\(.line // .original_line) | \(.body)"'
```

## Contents

- Figures
- Real examples with real parameters
- Documentation
- Numbers in tests
- Scope of the PR
- Structure and duplication
- No backward compatibility for unreleased formats
- CLI defaults
- Domain rules
- Efficiency

## Figures

She does not trust a test or a claim about sequences that is shown only as text.

| PR | Her words |
|---|---|
| #69 `evalue.rs` | "I'd really love figures for all these tests. It's hard to trust just pure text" |
| #67 `test_cli.rs` | "Why doesn't this grow to the left? Why only on the right? I need to see a figure" |
| #62 `search.rs` | "I want to see this example, show me the figure" |
| #58 `visualize_pair_html.py` | "I want to see an example output HTML too" |
| #67 | "Rename to docs/images/xdrop_walk_bcl2_ced9_bh1.png" and "reference this figure in some doc somewhere. it shouldn't be some bare file" |
| #73 `plot_ka_survival.py` | "Point people to this explainer ... This figure should be understandable to someone without context" |
| #54 `karlin_altschul.rs` | "Reference figures in the docs/images here" |

What to check:

- Each new or changed test about a region, an extension, a score or a fit has a figure of
  the case it tests, posted on the PR or committed under `docs/images/`.
- The file name says what and which example: `<what>_<protein pair>_<detail>.png`.
- Something links to the figure: a doc page, the README, or the doc comment of the code it
  explains. An unlinked image file is a finding.
- Read cold, the figure makes sense: title, axis labels, legend, and a link to the
  explainer page when it shows E-value or Karlin-Altschul quantities.
- A script that writes HTML or plots has an example of its output on the PR.
- Run the `clear-figures` checklist on every figure.

## Real examples with real parameters

| PR | Her words |
|---|---|
| #50 `test_cli.rs` | "I also want to see that the full BCL2/Ced9 region is still detected", with the region drawn as protein residues, HP letters and match bars |
| #54 `test_cli.rs` | "Show me the new BCL2/Ced9 region in this test" |
| #36 `index.rs` | "a test that actually shows the top k-mers by protein (k=5), dayhoff (k=7), and hp pbotc v1 (k=10) from the full CED9 protein sequence" |
| #62 `search.rs` | "k=9 for hp space ... is completely unreasonable. Show me at least k=12, but k=15+ is much more realistic. We usually start with k=19+ for hp ksizes" |
| #38 `search.rs` | "I want to see examples where the whole protein passes but the region doesn't" and "Use human-mouse gencode genes as an example" |
| #57 | "how does the expansion affect the pvalue, region enrichment, etc?" |
| #43, #56 | "Test ALL encodings here, not just a few", "Shouldn't this test all combos and not just a sliding pair?", "add tests of all other alphabets here" |
| #54 `karlin_altschul.rs` | "Why are you defining this inline? Don't we have a BCL2 variable you can use?" |
| #58 `pair_model.py` | "Always show the residues that are determined to be hydrophobic and polar for each alphabet, since it changes." |

The BCL2/Ced9 drawing she wants, from #50:

```
Ced9 pr: …RTVGNAQTD**QCPMSYGRLIGLISFGGFV**AAKMMESVE…
Ced9 hp: …pphhphppp**pphhphhphhhhhphhhhh**hhphhpphp…
                   |||||||||||||||||||
BCL2 hp: …hhphhpphh**pphhphhphhhhhphhhhh**phpphppph…
BCL2 pr: …FATVVEELF**RDGVNWGRIVAFFEFGGVM**CVESVNREM…
```

## Documentation

| PR | Her words |
|---|---|
| #54 `karlin_altschul.rs` | "Make sure ALL properties have a doc comment. Nothing should be empty" |
| #54 | "Suggest defaults for mismatch_penalty (2) and xdrop (8). Add doc comments" |
| #54 | "I feel like we should add some docs in general" |
| #54 `search.rs` | "wtf is r and wtf is t" |
| #68 `search.rs` | "what is hi?" and "Rename from a -> u as per the explainer" |
| #54 `test_cli.rs` | "Clarify what 'censored' means" |
| #38 `search.rs` | "what does non-modal mean?", "Explain multiplicity clearly and succinctly", "I don't understand what this is testing" |
| #71 `evalue.rs` | "Explain again what a decoy is? How is it defined, what does it represent" |
| #71 `main.rs` | "Explain how the closed form is defined", "I still don't understand --ka-null vs --ka-reference" |
| #38 `search.rs` | "Specify that lambda is the sum of database abundances of the k-mers from the region, divided by the region size, in the code and the PR text" |
| #54 `evalue.rs` | "State that the slope is K and intercept is lambda" |
| #38 | "this isn't a real pvalue ... So can we call it something else?" |
| #40 `index.rs` | "Add a comment explaining the reason for this", "I don't like this variable name, 'wanted'. Unclear WHY this is 'wanted'" |
| #38 `search.rs` | "what does saturating_sub help with here?" |
| #75 `docs/evalue.md` | "Is this file viewable on seanome.github.io/kmerseek?" |
| #50 | "Update changelog too" |

## Numbers in tests

| PR | Her words |
|---|---|
| #54 `karlin_altschul.rs` | "I really need to see the original citation to trust it. Give the originals and no magic numbers!" and "This applies to ALL the tests here" |
| #54 `test_cli.rs` | "Can the 242 be a variable instead of a magic number that we can use across multiple tests? Look for other re-used magic numbers" |
| #54 `index.rs` | "This feels like cheating. Why are you just putting these values in here?" |
| #38 `test_visualize_hits.py` | "Why are all these number fields strings? Aren't they integers and floats?" |

Tests still assert exact literals (see the `tests-assert-exact-values` memory). Each literal
also says where it came from: a cited paper and table, a hand calculation written out in the
comment, or the command that produced it.

## Scope of the PR

| PR | Her words |
|---|---|
| #54 `index.rs` | "Why are all these RocksDB changes in this PR?" and "Maybe that's a separate PR" |
| #38 `search.rs` | "Don't do this, removing low complexity in a separate PR" |

## Structure and duplication

| PR | Her words |
|---|---|
| #54 `main.rs` | "should this really be in main.rs? It's not the CLI" and "Move to karlin_altschul.rs" |
| #54 `search.rs` | "Why are there so many repeated properties? can't they be inherited or derived? This is hard to maintain" |
| #68 `search.rs` | "This seems stupid to have a single value in a struct? Shouldn't like r_database and m and n be here too?" |
| #53 `index.rs` | "Do we really need `pub fn load` in addition to `pub fn get_index_parameters`?" |
| #63 `index.rs` | "why are we duplicating things?" |
| #54 `evalue.rs`, #58 | "Do we really have to reinvent this wheel? Can't we just use a library?" then, for biopython, "weigh the pros and cons of having such a big dependency" |
| #37 `index.rs` | "Is this necessary? Can't we just make the property public?" |

## No backward compatibility for unreleased formats

| PR | Her words |
|---|---|
| #53 `index.rs` | "I don't care about the old indices. This package doesn't have a ton of users. I care about making the future better" |
| #43 | "Don't reference older indices", "Don't reference 'older spellings' let's make this a clean break", "nothing is 'old' here" |
| #37 `index.rs` | "Can we add the kmerseek version to the index and read the field that way?" |

The one exception she asked for: names sourmash uses (`hp`, `protein`, `dayhoff`) keep
parsing, for compatibility with sourmash (#43).

## CLI defaults

| PR | Her words |
|---|---|
| #68 `main.rs` | "Why is this required? It should just be used directly ... by reading the index ... which should be the default" |
| #71 `main.rs` | "does running `kmerseek index` or `kmerseek search` require `--ka-k`? Why can't it use the value from the index?" and "Why doesn't running `kmerseek index` produce these fits?" |
| #70 `evalue.rs` | "Can we default to ShuffledDipeptide?" |
| #37 `main.rs` | "Let's also add --remove-low-complexity to the search tool" |

## Domain rules

| PR | Her words |
|---|---|
| #63 `search.rs` | "always skip self-matches ... increases the row count by N" |
| #74 `search.rs` | "do you deal with the ambiguous k-mers properly, using total windows and not total k-mers?" |
| #49 `aminoacid.rs` | cap ambiguous residues per window; she asked for "the densest window in UniRef50/90/100" before choosing the cap |

## Efficiency

| PR | Her words |
|---|---|
| #37 `index.rs` | "is this necessary to re-iterate over all the signatures after building them? Couldn't we return this info as they're built?" |
| #36 `index.rs` | "could we just store the top N and bottom N k-mers as we go along? I'd rather not add any more O(N) or bigger operations" |
| #36 `index.rs` | "Why do we need to do resolve_kmer_string? Can't we use the saved kmers from the index?" |
| #48 `index.rs` | "This feels like a race condition waiting to happen. Let's avoid" |
| #37 `main.rs` | "how much time does removing low complexity add to index and search?" |
| #38 `search.rs` | "Is this the most efficient iteration? What kind of time are we losing here?" |
