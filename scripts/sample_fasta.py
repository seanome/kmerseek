#!/usr/bin/env python3
"""Draw a seeded random sample of sequences from a FASTA, for the Karlin-Altschul fits
in docs/evalue.md.

    python scripts/sample_fasta.py uniprot_sprot.fasta.gz 15000 --seed 1 -o swissprot_sample15000.fasta

Reads the file twice: once to count the records, once to write the ones picked. The
records picked are random.Random(seed).sample(range(count), n), written in file order,
so the same file, n and seed always give the same sample. Prints the count, the residues
written and the input's MD5, which docs/evalue.md records next to the fit.
"""

import argparse
import gzip
import hashlib
import random


def records(path):
    """(header line, sequence) pairs, in file order."""
    opener = gzip.open if path.endswith(".gz") else open
    header, seq = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line, []
            else:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def md5(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("fasta")
    parser.add_argument("n", type=int)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    count = sum(1 for _ in records(args.fasta))
    picked = set(random.Random(args.seed).sample(range(count), args.n))
    residues = 0
    with open(args.output, "w") as out:
        for i, (header, seq) in enumerate(records(args.fasta)):
            if i in picked:
                out.write(f"{header}{seq}\n")
                residues += len(seq)
    print(f"{args.fasta}\tmd5 {md5(args.fasta)}\t{count} records\t"
          f"{args.n} sampled with seed {args.seed}\t{residues} residues\t-> {args.output}")


if __name__ == "__main__":
    main()
