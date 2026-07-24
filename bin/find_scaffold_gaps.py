#!/usr/bin/env python3

import argparse
import re

parser = \
    argparse.ArgumentParser(

        description = "Returns positions of gaps in a \
            scaffolded genome assembly."

    )

parser.add_argument(

    "fasta",
    type = str,
    help = "Path to genome assembly in FASTA format."

)

options = parser.parse_args()

def parse_fasta(fasta):

    fasta_dict = {}

    with open(fasta, "r") as f:

        header = "" 
 
        for line in f:

            if line.startswith(">"):
                header = line.strip()[1:]
                fasta_dict[header] = ""

            else:
                
                if header:
                    fasta_dict[header] += line.strip()

                else:
                    raise ValueError(
                        "FASTA file has sequence before header."
                    )
        
    return fasta_dict

def find_gaps(fasta_dict):

    for header, seq in fasta_dict.items():

        for match in re.finditer(r'[Nn]+', seq):
            start = match.start()
            end = match.end()
            print(f"{header}\t{start}\t{end}")

if __name__ == "__main__":
    fasta_dict = parse_fasta(options.fasta)
    find_gaps(fasta_dict)