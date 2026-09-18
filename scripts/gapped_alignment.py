"""Gapped alignment of two short protein stretches, for counting identities.

A run of shared k-mers is an exact match on one diagonal in the reduced alphabet. Two
proteins that are homologous there often align a few residues off that diagonal once
gaps are allowed (MCL-1's BH1 run reads 1 identical residue on its seed diagonal though
NWGR is in both rows), so identities are counted after a gapped alignment of the
run and its flanks: BLOSUM62, gap open 11, gap extend 1, the BLAST defaults.

`local_align(a, b)` returns the best local (Smith-Waterman) alignment; `global_align(a, b)`
the best end-to-end (Needleman-Wunsch) alignment, which is what a run and its flanks get,
since a local alignment of same-class but different residues shrinks to its few identical
ones. Both return 0-based half-open coordinates on `a` and `b` and the two aligned rows
with `-` for gaps. Pure numpy: the row loop is over `a`, and the horizontal-gap state
along `b` is a running maximum, so nothing is quadratic in Python.
"""

import numpy as np

_AA = "ARNDCQEGHILKMFPSTWYVBZX*"
_BLOSUM62_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""
BLOSUM62 = np.array([[int(v) for v in line.split()] for line in _BLOSUM62_ROWS.strip().splitlines()], dtype=np.int32)
_INDEX = {aa: i for i, aa in enumerate(_AA)}
GAP_OPEN, GAP_EXTEND = 11, 1
NEG = -(10**6)


def _codes(seq):
    """Residue indices into BLOSUM62; anything unknown (U, O, lower case) scores as X."""
    return np.array([_INDEX.get(c.upper(), _INDEX["X"]) for c in seq], dtype=np.intp)


def _fill(a, b, local):
    """Gotoh matrices: M (match ends here), X (gap in `a`, i.e. a run of `b` residues against
    `-`), Y (gap in `b`). Row i is `a[i-1]`, column j is `b[j-1]`. `local` clamps M at zero
    and makes leading gaps free; global mode charges them like any other gap."""
    n, m = len(a), len(b)
    ca, cb = _codes(a), _codes(b)
    M = np.full((n + 1, m + 1), 0 if local else NEG, dtype=np.int32)
    X = np.full((n + 1, m + 1), NEG, dtype=np.int32)
    Y = np.full((n + 1, m + 1), NEG, dtype=np.int32)
    js = np.arange(1, m + 1, dtype=np.int32)
    M[0, 0] = 0
    if not local:
        X[0, 1:] = -GAP_OPEN - GAP_EXTEND * (js - 1)
        Y[1:, 0] = -GAP_OPEN - GAP_EXTEND * (np.arange(n, dtype=np.int32))
    for i in range(1, n + 1):
        sub = BLOSUM62[ca[i - 1], cb]
        prev_best = np.maximum(np.maximum(M[i - 1, :-1], X[i - 1, :-1]), Y[i - 1, :-1])
        M[i, 1:] = np.maximum(0, prev_best + sub) if local else prev_best + sub
        # A vertical gap (a residue of `a` against `-`) extends from the row above.
        Y[i, 1:] = np.maximum(M[i - 1, 1:] - GAP_OPEN, Y[i - 1, 1:] - GAP_EXTEND)
        # A horizontal gap ends at column j after opening at some k < j: a running maximum
        # of M[i, k] - open + extend * (k + 1) minus extend * j, with X's own extension folded
        # in because the max already carries the cheapest way to reach each column. In global
        # mode a gap may also open from Y (a gap in the other sequence just closed).
        opener = M[i, :-1] if local else np.maximum(M[i, :-1], Y[i, :-1])
        opened = opener - GAP_OPEN + GAP_EXTEND * js
        X[i, 1:] = np.maximum.accumulate(opened) - GAP_EXTEND * js
        if not local:
            Y[i, 1:] = np.maximum(Y[i, 1:], X[i - 1, 1:] - GAP_OPEN)
    return M, X, Y


def _traceback(a, b, M, X, Y, i, j, state, local):
    """Walk back from cell (i, j) in `state`, returning the aligned rows and the start cell."""
    ra, rb = [], []
    while i > 0 or j > 0:
        if i == 0:
            state = "X"
        elif j == 0:
            state = "Y"
        if state == "M":
            if local and M[i, j] == 0:
                break
            ra.append(a[i - 1])
            rb.append(b[j - 1])
            sub = BLOSUM62[_INDEX.get(a[i - 1].upper(), _INDEX["X"]), _INDEX.get(b[j - 1].upper(), _INDEX["X"])]
            prev = M[i, j] - sub
            i, j = i - 1, j - 1
            state = "M" if M[i, j] == prev else ("X" if X[i, j] == prev else "Y")
        elif state == "Y":
            ra.append(a[i - 1])
            rb.append("-")
            came_from_m = M[i - 1, j] - GAP_OPEN == Y[i, j]
            came_from_x = not local and X[i - 1, j] - GAP_OPEN == Y[i, j]
            state = "M" if came_from_m else ("X" if came_from_x else "Y")
            i -= 1
        else:
            ra.append("-")
            rb.append(b[j - 1])
            came_from_m = M[i, j - 1] - GAP_OPEN == X[i, j]
            came_from_y = not local and Y[i, j - 1] - GAP_OPEN == X[i, j]
            state = "M" if came_from_m else ("Y" if came_from_y else "X")
            j -= 1
        if local and i == 0 or local and j == 0:
            break
    return "".join(reversed(ra)), "".join(reversed(rb)), i, j


def _result(a, b, M, X, Y, end_i, end_j, state, local):
    row_a, row_b, start_i, start_j = _traceback(a, b, M, X, Y, int(end_i), int(end_j), state, local)
    score = {"M": M, "X": X, "Y": Y}[state][end_i, end_j]
    return {
        "score": int(score),
        "a_start": start_i,
        "a_end": int(end_i),
        "b_start": start_j,
        "b_end": int(end_j),
        "a_row": row_a,
        "b_row": row_b,
        "identical": sum(x == y and x != "-" for x, y in zip(row_a, row_b)),
        "aligned": sum(x != "-" and y != "-" for x, y in zip(row_a, row_b)),
        "columns": len(row_a),
    }


def local_align(a, b):
    """Best Smith-Waterman alignment of `a` and `b`. Returns None when nothing scores above
    zero. Coordinates are 0-based half-open on the inputs."""
    if not a or not b:
        return None
    M, X, Y = _fill(a, b, local=True)
    end_i, end_j = np.unravel_index(np.argmax(M), M.shape)
    if M[end_i, end_j] <= 0:
        return None
    return _result(a, b, M, X, Y, end_i, end_j, "M", local=True)


def global_align(a, b):
    """Best end-to-end alignment of `a` and `b`, end gaps charged like any other. Returns
    None for an empty input."""
    if not a or not b:
        return None
    M, X, Y = _fill(a, b, local=False)
    n, m = len(a), len(b)
    state = max(("M", "X", "Y"), key=lambda st: {"M": M, "X": X, "Y": Y}[st][n, m])
    return _result(a, b, M, X, Y, n, m, state, local=False)
