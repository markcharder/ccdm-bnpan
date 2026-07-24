#!/usr/bin/env python3

import argparse
import sys

parser = argparse.ArgumentParser(
    description = "Remove contigs below a certain size from a FASTA."
)

parser.add_argument(
    "--fasta",
    "-f",
    required=True,
    help="Input FASTA file."
)

parser.add_argument(
    "--min-size",
    "-m",
    type=int,
    default=1000,
    help="Minimum contig size to retain (default: 1000)."
)

def filter_fasta(
        
    fasta_file: str,
    min_size: int
    ):
    fasta_seqs = dict()
    header = ""
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()
                fasta_seqs[header] = ""
            else:
                if header:
                    fasta_seqs[header] += line.strip()
                else:
                    sys.stderr.write("Malformed FASTA. Sequence line with no header.")
                    sys.exit(1)

    for header, seq in fasta_seqs.items():
        if len(seq) >= min_size:
            print(header)
            print(seq)

if __name__ == "__main__":
    args = parser.parse_args()
    filter_fasta(
        fasta_file=args.fasta,
        min_size=args.min_size
    )