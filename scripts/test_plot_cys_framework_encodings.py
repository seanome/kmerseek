"""Pins the partitions in plot_cys_framework_encodings.py to src/rust/alphabets.rs.

The script writes the h/p/c residue classes out as strings. These tests read the
`build_hp` / `build_hpc` calls in the Rust file and fail if the strings drift, so the
figure and the README counts cannot go stale without a test failing.

Run with: python -m pytest scripts/test_plot_cys_framework_encodings.py -v
"""

import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
import plot_cys_framework_encodings as fig

ALPHABETS_RS = os.path.join(os.path.dirname(__file__), "..", "src", "rust", "alphabets.rs")


def rust_partition(static_name):
    """The byte-string arguments of the build_hp/build_hpc call that defines `static_name`."""
    with open(ALPHABETS_RS) as handle:
        source = handle.read()
    match = re.search(
        rf"static {static_name}\b.*?build_hpc?\(([^)]*)\)", source, flags=re.DOTALL
    )
    assert match, f"no build_hp/build_hpc call found for {static_name}"
    return re.findall(r'b"([A-Z]+)"', match.group(1))


def test_pbotc_two_class_partition_matches_rust():
    h, p = rust_partition("PBOTC_1ST_ED_HP")
    assert (fig.PBOTC_H, fig.PBOTC_P) == (h, p)


def test_pbotc_three_class_partition_matches_rust():
    h, p, c = rust_partition("PBOTC_1ST_ED_HPC")
    assert (fig.HPC_H, fig.HPC_P, "C") == (h, p, c)
