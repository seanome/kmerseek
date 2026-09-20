"""Gapped alignment of two short protein stretches, for counting identities.

A run of shared k-mers is an exact match on one diagonal in the reduced alphabet. Two
proteins that are homologous there often align a few residues off that diagonal once
gaps are allowed (MCL-1's BH1 run reads 1 identical residue on its seed diagonal though
NWGR is in both rows), so identities are counted after a gapped alignment of the run and
its flanks. Biopython's PairwiseAligner does the alignment with BLAST's protein defaults:
BLOSUM62, gap open 11, gap extend 1.
"""

from Bio.Align import PairwiseAligner, substitution_matrices

_BLOSUM62 = substitution_matrices.load("BLOSUM62")
_ALIGNER = PairwiseAligner(substitution_matrix=_BLOSUM62, open_gap_score=-11, extend_gap_score=-1, mode="global")


def _scorable(seq):
    """Upper case, with any residue BLOSUM62 has no row for (U, O) scored as X."""
    return "".join(c if c in _BLOSUM62.alphabet else "X" for c in seq.upper())


def _row(aligned, original):
    """The aligned row with the caller's own letters back in place of the scorable ones."""
    letters = iter(original)
    return "".join("-" if c == "-" else next(letters) for c in aligned)


def global_align(a, b):
    """Best end-to-end alignment of `a` and `b`, end gaps charged like any other. Returns
    None for an empty input. Coordinates are 0-based half-open on the inputs; the rows
    carry `-` for gaps."""
    if not a or not b:
        return None
    hit = _ALIGNER.align(_scorable(a), _scorable(b))[0]
    row_a, row_b = _row(hit[0], a), _row(hit[1], b)
    (a_start, a_end), (b_start, b_end) = (hit.coordinates[i][[0, -1]] for i in (0, 1))
    return {
        "score": int(hit.score),
        "a_start": int(a_start),
        "a_end": int(a_end),
        "b_start": int(b_start),
        "b_end": int(b_end),
        "a_row": row_a,
        "b_row": row_b,
        "identical": sum(x == y and x != "-" for x, y in zip(row_a, row_b)),
        "aligned": sum(x != "-" and y != "-" for x, y in zip(row_a, row_b)),
        "columns": len(row_a),
    }
