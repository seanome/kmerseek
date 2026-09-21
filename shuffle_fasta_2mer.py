#!/usr/bin/env python3
"""Shuffle FASTA sequences while keeping every 2-mer count, so a shuffled protein has the
same dipeptide composition as the original but no other resemblance to it. Writes N
shuffled versions of each sequence, named `<header>_shuffle<i>`.

Uses ushuffle when it is installed and otherwise the same algorithm in Python (Altschul
and Erickson 1985, the one behind ushuffle): pick, for every residue but the last, the
transition it will leave on last, keep picking until those last transitions all lead to
the final residue, then walk the transitions in a random order.

    python shuffle_fasta_2mer.py targets.fasta -n 20 --seed 1 -o decoys.fasta
"""

import argparse
import random
import sys
from collections import defaultdict

try:
    import ushuffle
except ImportError:
    ushuffle = None


def parse_fasta(filepath):
    """Simple FASTA parser yielding (header, sequence) tuples."""
    header = None
    seq_lines = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_lines)
                header = line[1:]  # Remove '>'
                seq_lines = []
            else:
                seq_lines.append(line)

        # Don't forget the last sequence
        if header is not None:
            yield header, ''.join(seq_lines)


def _last_edges_reach_end(last_edge, end):
    """Whether following each residue's chosen last transition always arrives at `end`."""
    for start in last_edge:
        seen, v = set(), start
        while v != end:
            if v in seen:
                return False
            seen.add(v)
            v = last_edge[v]
    return True


def shuffle_2mer_python(seq, rng):
    """A random sequence with exactly the 2-mer counts of `seq`."""
    if len(seq) < 3:
        return seq
    edges = defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        edges[a].append(b)
    end = seq[-1]
    while True:
        last_edge = {v: rng.choice(out) for v, out in edges.items() if v != end}
        if _last_edges_reach_end(last_edge, end):
            break
    order = {}
    for v, out in edges.items():
        rest = list(out)
        if v != end:
            rest.remove(last_edge[v])
        rng.shuffle(rest)
        order[v] = rest + ([last_edge[v]] if v != end else [])
    out, v = [seq[0]], seq[0]
    for _ in range(len(seq) - 1):
        v = order[v].pop(0)
        out.append(v)
    return ''.join(out)


def shuffle_sequence(seq, k=2, rng=None):
    """Shuffle sequence preserving k-mer composition. The Python fallback does k=2 only."""
    if ushuffle is not None:
        return ushuffle.shuffle(seq, len(seq), k)
    if k != 2:
        sys.exit(f"k={k} needs ushuffle (pip install ushuffle); the built-in shuffle keeps 2-mers only")
    return shuffle_2mer_python(seq, rng or random.Random())


def main():
    parser = argparse.ArgumentParser(
        description='Shuffle FASTA sequences preserving 2-mer composition'
    )
    parser.add_argument('input', help='Input FASTA file')
    parser.add_argument('-o', '--output', default='-',
                        help='Output FASTA file (default: stdout)')
    parser.add_argument('-n', '--num-shuffles', type=int, default=10,
                        help='Number of shuffled versions per sequence (default: 10)')
    parser.add_argument('-k', '--kmer-size', type=int, default=2,
                        help='K-mer size to preserve (default: 2)')
    parser.add_argument('--seed', type=int,
                        help='random seed, for a reproducible set of shuffles')
    args = parser.parse_args()

    if args.seed is not None:
        if ushuffle is not None:
            ushuffle.set_seed(args.seed)
    rng = random.Random(args.seed)
    out_fh = sys.stdout if args.output == '-' else open(args.output, 'w')

    try:
        for header, seq in parse_fasta(args.input):
            for i in range(1, args.num_shuffles + 1):
                shuffled_seq = shuffle_sequence(seq, k=args.kmer_size, rng=rng)
                # Write shuffled sequence with modified header
                out_fh.write(f'>{header}_shuffle{i}\n')
                # Write sequence in 60-character lines
                for j in range(0, len(shuffled_seq), 60):
                    out_fh.write(shuffled_seq[j:j+60] + '\n')
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()


if __name__ == '__main__':
    main()
