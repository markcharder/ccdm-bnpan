#!/usr/bin/env python3

import argparse
"""
Takes the VCF output from vg deconstruct and does this:
1) Sets genotypes to missing for samples in the list of conflicts.
2) Removes the CONFLICT field from the INFO column and header.
3) Converts haploid genotypes to diploid by duplicating the allele.
4) Adds samples to the header.
"""
def parse_conflicts(
        info: str,
        samples: list,
        genotypes: list
        ):

    conflicting_samples = \
        set(info.split("CONFLICT=")[1].split(";")[0].split(","))

    for i, sample in enumerate(samples):

        if sample in conflicting_samples:
            genotypes[i] = "./."
        
        else:
            genotypes[i] = genotypes[i]+"/"+genotypes[i]

    return(genotypes)

def fix_conflicts(vcf: str):

    with open(vcf, 'r') as f:
        samples = ""

        for line in f:

            if line.startswith("##"):
                
                if "CONFLICT" in line:
                    continue

                else:
                    print(line.strip())
                    continue

            if line.startswith("#CHROM"):
                parts = line.strip().split("\t")
                samples = parts[9:]
                samples_header = \
                    "\n".join(["##SAMPLE=<ID={0}>".format(s) for s in samples])
                print(samples_header)
                print(line.strip())
                continue

            parts = line.strip().split("\t")
            info = parts[7]

            if samples:

                if "CONFLICT=" not in info:
                    genotypes = [ gt + "/" + gt for gt in parts[9:] ]
                    parts[9:] = genotypes
                    print("\t".join(parts))
                    continue
                else:
                    genotypes = parts[9:]
                    genotypes = \
                        parse_conflicts(info, samples, genotypes)                
                    parts[7] = info.split(";CONFLICT=")[0]
                    parts[9:] = genotypes
                    print("\t".join(parts))

            else:
                exit("Error: No samples found in VCF.")
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="VCF for fixing.")
    options = parser.parse_args()
    fix_conflicts(options.vcf)