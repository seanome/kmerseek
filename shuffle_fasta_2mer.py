#!/usr/bin/env python3
"""
Shuffle FASTA sequences while preserving 2-mer composition using ushuffle.
Generates N shuffled versions of each sequence.
"""

import argparse
import sys
import ushuffle


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


def shuffle_sequence(seq, k=2):
    """Shuffle sequence preserving k-mer composition."""
    shuffled = ushuffle.shuffle(seq, len(seq), k)
    return shuffled


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
    args = parser.parse_args()

    out_fh = sys.stdout if args.output == '-' else open(args.output, 'w')

    try:
        for header, seq in parse_fasta(args.input):
            for i in range(1, args.num_shuffles + 1):
                shuffled_seq = shuffle_sequence(seq, k=args.kmer_size)
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
