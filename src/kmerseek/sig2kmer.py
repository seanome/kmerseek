# ! /usr/bin/env python3
"""
Given a signature file and a collection of sequences, output all of the
k-mers and sequences that match a hashval in the signature file.

Cribbed from https://github.com/dib-lab/sourmash/pull/724/
"""

import sourmash
from sourmash import MinHash
from sourmash.minhash import hash_murmur
import screed
import csv
from sourmash.logging import notify


NOTIFY_EVERY_BP = int(1e5)


def degenerate_protein_chatgpt(sequence, moltype):
    """Convert protein sequence from 20-letter amino acid alphabet to degenerate alphabet"""
    # Pre-determine the alphabet function based on moltype
    if moltype == "hp":
        alphabet_func = sourmash._lowlevel.lib.sourmash_aa_to_hp
    elif moltype == "dayhoff":
        alphabet_func = sourmash._lowlevel.lib.sourmash_aa_to_dayhoff
    else:
        raise ValueError(f"Unknown moltype: {moltype}")

    # Convert the entire sequence to bytes once
    byte_sequence = sequence.encode("utf-8")

    # Apply the alphabet function to each byte in the sequence and join the result
    degenerate = b"".join(
        alphabet_func(letter.to_bytes(1, "big")) for letter in byte_sequence
    ).decode()

    return degenerate


def degenerate_protein(sequence, moltype):
    """Convert protein sequence from 20-letter amino acid alphabet to degenerate alphabet"""
    if moltype == "hp":
        alphabet = sourmash._lowlevel.lib.sourmash_aa_to_hp
    elif moltype == "dayhoff":
        alphabet = sourmash._lowlevel.lib.sourmash_aa_to_dayhoff
    else:
        raise ValueError(f"Unknown moltype: {moltype}")

    # Convert the entire sequence to bytes once
    byte_encoded = (x.encode("utf-8") for x in sequence)

    degenerate = b"".join(alphabet(letter) for letter in byte_encoded).decode()
    return degenerate


def get_kmer_moltype(sequence, start, ksize, moltype, input_is_protein):
    kmer_in_seq = sequence[start : start + ksize]
    if moltype == "DNA":
        # Get reverse complement
        kmer_rc = screed.rc(kmer_in_seq)
        if kmer_in_seq > kmer_rc:  # choose fwd or rc
            kmer_encoded = kmer_rc
        else:
            kmer_encoded = kmer_in_seq
    elif input_is_protein:
        kmer_encoded = degenerate_protein(kmer_in_seq, moltype)
    elif not input_is_protein:
        raise NotImplementedError("Currently cannot translate DNA to protein sequence")
    return kmer_encoded, kmer_in_seq


def revise_ksize(ksize, moltype, input_is_protein):
    """If input is protein, then divide the ksize by three"""
    if moltype == "DNA":
        return ksize
    elif input_is_protein:
        # Ksize includes codons
        return int(ksize / 3)
    else:
        return ksize


def get_kmers_for_hashvals(
    sequence, hashvals, ksize, moltype, input_is_protein, n_kmers_in_sequence, minhash
):
    "Return k-mers from 'sequence' that yield hashes in 'hashvals'."
    # uppercase!
    sequence = sequence.upper()

    # Divide ksize by 3 if sequence is protein
    # ksize = revise_ksize(ksize, moltype, input_is_protein)

    for start in range(0, len(sequence) - ksize + 1):

        kmer_encoded, kmer_in_seq = get_kmer_moltype(
            sequence, start, ksize, moltype, input_is_protein
        )

        # NOTE: we do not avoid non-ACGT characters, because those k-mers,
        # when hashed, shouldn't match anything that sourmash outputs.
        hashval = hash_murmur(kmer_encoded)
        print(
            f"kmer_encoded: {kmer_encoded}, kmer_in_seq: {kmer_in_seq}, hashval: {hashval}"
        )
        import pdb

        pdb.set_trace()
        if hashval in hashvals:
            n_kmers_in_sequence += 1
            yield kmer_encoded, kmer_in_seq, hashval, start


def get_matching_hashes_in_file(
    filename,
    ksize,
    moltype,
    minhash,
    input_is_protein,
    hashes,
    found_kmers,
    m,
    n,
    n_seq,
    seqout_fp,
    kmerout_w,
    watermark,
    first=False,
):
    for record in screed.open(filename):
        n += len(record.sequence)
        n_seq += 1
        while n >= watermark:
            notify(
                "...Searched {:d} residues,\tfound {} kmers in\t{} seqs from\t{}",
                watermark,
                found_kmers,
                n_seq,
                filename,
                end="\r",
            )
            watermark += NOTIFY_EVERY_BP

        n_kmers_in_sequence = 0
        # now do the hard work of finding the matching k-mers!
        for kmer_encoded, kmer_in_seq, hashval, i in get_kmers_for_hashvals(
            record.sequence,
            hashes,
            ksize,
            moltype,
            input_is_protein,
            n_kmers_in_sequence,
            minhash,
        ):
            found_kmers += 1

            if n_kmers_in_sequence != len(hashes):
                raise ValueError(
                    f"Did not find all hashvals in sequence k-mers. Found "
                    f"only {n_kmers_in_sequence} k-mers out of {len(hashes)}"
                )

            # write out sequence
            if seqout_fp:
                seqout_fp.write(
                    ">{}|hashval:{}|kmer:{}|kmer_encoded:{}\n{}\n".format(
                        record.name, hashval, kmer_in_seq, kmer_encoded, record.sequence
                    )
                )
                m += len(record.sequence)
            if kmerout_w:
                kmerout_w.writerow(
                    [
                        kmer_in_seq,
                        kmer_encoded,
                        str(hashval),
                        i,
                        record["name"],
                        filename,
                    ]
                )

            n_kmers_in_sequence = 0
            if first:
                return m, n
    return m, n


def _make_kmer_filename(sig):
    return f"{sig}.kmers.csv"


def _setup_kmer_writer(sig):
    output_kmers = _make_kmer_filename(sig)
    kmerout_fp = open(output_kmers, "w")
    kmerout_w = csv.writer(kmerout_fp)
    kmerout_w.writerow(
        [
            "kmer_in_sequence",
            "kmer_in_alphabet",
            "hashval",
            "start",
            "read_name",
            "filename",
        ]
    )
    return output_kmers, kmerout_fp, kmerout_w


def get_kmers(sig, fasta, moltype, ksize, scaled):
    # Scale is not used, but is accepted for ease of keyword argument sharing
    print("\n--- get_kmers ---")
    print(f"Getting all k-mers from file {sig}")

    output_kmers, kmerout_fp, kmerout_w = _setup_kmer_writer(sig)

    # first, load the signature and extract the hashvals
    sigobj = list(
        sourmash.load_file_as_signatures(sig, ksize=ksize, select_moltype=moltype)
    ).pop()
    query_hashvals = set(sigobj.minhash.hashes.keys())
    query_ksize = sigobj.minhash.ksize
    new_minhash = MinHash(
        n=0,
        ksize=sigobj.minhash.ksize,
        hp=sigobj.minhash.hp,
        dayhoff=sigobj.minhash.dayhoff,
        is_protein=sigobj.minhash.is_protein,
        scaled=scaled,
    )

    import pdb

    pdb.set_trace()

    print(f"query_ksize: {query_ksize}")
    print(f"ksize: {ksize}")

    # TODO: Simplify and remove all the tracking numbers, e.g. number residues loaded, kmers found, etc.
    # I really don't like how this relies on so many parameters for tracking how far we've
    # progressed through the files. I wish this function was massively simpler, not accepting any of the
    # keyword arguments initialized with a 0
    n = 0
    found_kmers = 0
    n_seq = 0
    m = 0
    get_matching_hashes_in_file(
        fasta,
        query_ksize,
        moltype,
        new_minhash,
        input_is_protein=True,
        hashes=query_hashvals,
        found_kmers=0,
        m=m,  # bp in found sequences
        n=n,  # bp loaded
        n_seq=n_seq,
        seqout_fp=None,
        kmerout_w=kmerout_w,
        watermark=int(NOTIFY_EVERY_BP),
    )

    notify(
        "read {} bp, found {} kmers matching hashvals\nwrote to {}",
        n,
        found_kmers,
        output_kmers,
    )
    return output_kmers


class Args:
    """Dummy class to call what is normally a command line function from Python"""

    def __init__(self, sig, fasta, moltype, ksize, scaled):
        self.signatures = [sig]  # List of signature files
        self.sequences = [fasta]  # List of sequence files
        self.ksize = ksize
        self.quiet = False
        self.force = False
        self.moltype = moltype
        self.translate = False
        self.check_sequence = True
        self.save_kmers = None
        self.save_sequences = None
        self.picklist = None
        self.scaled = scaled
        self.dna = False
        self.skipm1n3 = False
        self.skipm2n3 = False
        self.from_file = False

        self.dayhoff = False
        self.hp = False
        self.protein = False
        if moltype == "protein":
            self.protein = True
        elif moltype == "dayhoff":
            self.dayhoff = True
        elif moltype == "hp":
            self.hp = True


def get_kmers_cli(sig, fasta, moltype, ksize, scaled):
    # Example usage:
    args = Args(
        sig=sig,
        fasta=fasta,
        moltype=moltype,  # or "dna", "dayhoff", "hp"
        ksize=ksize,
        scaled=scaled,
    )

    sourmash.sig.__main__.kmers(args)
