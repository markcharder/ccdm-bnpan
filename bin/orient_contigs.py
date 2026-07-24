#!/usr/bin/env python3
"""
Orient query contigs relative to a reference using a PAF alignment.

For each query contig, sums the aligned bases on + and - strand.
If more bases align on -, the contig is reverse-complemented in the output FASTA.

Usage:
    orient_contigs.py <query.fa> <alignments.paf> <output.fa>
"""

import sys
from collections import defaultdict


def reverse_complement(seq):
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def read_fasta(path):
    sequences = {}
    name = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        sequences[name] = "".join(buf)
    return sequences


def parse_paf(path):
    """Return per-contig aligned base counts by strand."""
    plus_bases  = defaultdict(int)
    minus_bases = defaultdict(int)
    with open(path) as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) < 12:
                continue
            qname  = parts[0]
            alnlen = int(parts[10])
            strand = parts[4]
            if strand == "+":
                plus_bases[qname]  += alnlen
            else:
                minus_bases[qname] += alnlen
    return plus_bases, minus_bases


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    query_fa, paf_file, out_fa = sys.argv[1], sys.argv[2], sys.argv[3]

    print("Reading query FASTA...", file=sys.stderr)
    sequences = read_fasta(query_fa)

    print("Parsing PAF...", file=sys.stderr)
    plus_bases, minus_bases = parse_paf(paf_file)

    to_rc = set()
    for name in sequences:
        p = plus_bases.get(name, 0)
        m = minus_bases.get(name, 0)
        if m > p:
            to_rc.add(name)
            print(f"  RC: {name}  (+ {p} bp  - {m} bp)", file=sys.stderr)
        else:
            print(f"  OK: {name}  (+ {p} bp  - {m} bp)", file=sys.stderr)

    print(f"\nReverse-complementing {len(to_rc)} of {len(sequences)} contigs",
          file=sys.stderr)

    with open(out_fa, "w") as fh:
        for name, seq in sequences.items():
            if name in to_rc:
                fh.write(f">{name}_rc\n")
                fh.write(reverse_complement(seq) + "\n")
            else:
                fh.write(f">{name}\n")
                fh.write(seq + "\n")

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
