#!/usr/bin/env python3

import argparse
import sys

def parse_families(families):

    fam_dict = dict()

    with open(families, 'r') as f:
        
        for line in f:
            subfam, fam = line.strip().split('\t')
            fam_dict[subfam] = fam

    return fam_dict

def annotate_gff(gff, fam_dict):

    with open(gff, 'r') as f:

        for line in f:
            parts = line.strip().split('\t')
            subfam = parts[8].split(":")[1].split(" ")[0][:-1]

            if subfam not in fam_dict:
                sys.stderr.write(subfam + "\tnot found\n")
                print(line.strip())
            
            else:
                line = line.replace(line.strip(), line.strip() + ";FAM=" + fam_dict[subfam])
                print(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--families",
        help = "Path to the families TSV file. Col1: sub-family name, Col2: family.",
        required = True
        )

    parser.add_argument(
        "-g", "--gff",
        help = "GFF file to annotate.",
        required = True
        )
    args = parser.parse_args()
    fam_dict = parse_families(args.families)
    annotate_gff(args.gff, fam_dict)
