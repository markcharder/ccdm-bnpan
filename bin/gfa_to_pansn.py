#!/usr/bin/env python3

import sys
import re

def convert_gfa(input_file, output_file):

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:

        for line in infile:

            if line.startswith('W\t'):
                parts = line.strip().split('\t')                
                sample = parts[1]
                haplotype = parts[2]
                contig = parts[3]                
                path_name = f"{sample}#{haplotype}#{contig}"
                walk_str = parts[6]
                steps = re.findall(r'([><])([^><]+)', walk_str)
                p_path = ",".join(
                        [
                        f"{node}{'+' if direction == '>' else '-'}"
                        for direction, node in steps
                        ]
                    )                
                overlaps = "*"
                outfile.write(f"P\t{path_name}\t{p_path}\t{overlaps}\n")

            else:
                outfile.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_w_to_p.py <input.gfa> <output.gfa>")
        sys.exit(1)
        
    convert_gfa(sys.argv[1], sys.argv[2])