#!/usr/bin/env python3

import csv
import glob
import os
import re

READS_DIR = "data/graphmap/reads_illumina"
CSV_PATH = "doc/info/sample-ids.csv"

def find_pair(prefix: str):
    if not prefix or prefix == "NONE" or prefix == "NA":
        return []
    pattern = os.path.join(READS_DIR, f"{prefix}_*")
    hits = glob.glob(pattern)
    return hits

def sort_pair(paths):

    def key(p):
        base = os.path.basename(p)
        if re.search(r"(_1\b|_1\.)", base): return 0
        if re.search(r"(_f1\b|_f1\.)", base): return 0
        if re.search(r"(_R1_trimmed\b|_R1_trimmed\.)", base): return 0

        if re.search(r"(_2\b|_2\.)", base): return 1
        if re.search(r"(_r2\b|_r2\.)", base): return 1
        if re.search(r"(_R2_trimmed\b|_R2_trimmed\.)", base): return 1

        return 2

    return sorted(paths, key=key)

with open(CSV_PATH, newline="") as fh:

    reader = csv.DictReader(fh)

    i=0

    with open("data/graphmap/conf/variety_ids.txt", 'w') as f1, open("data/graphmap/conf/variety_ids_conv.txt", 'w') as f2:

        for row in reader:
            i+=1
            variety = row.get("VARIETY", "")
            variety_code = "V" + str(i) + "SR"
            curtin = row.get("CURTIN_ID", "").strip()
            seqacc = row.get("SEQUENCING_ACCESSION", "").strip()

            matches = []

            for prefix in (seqacc, curtin):
                matches = find_pair(prefix)

                if matches:
                
                    if len(matches) == 2:
                        r1, r2 = sort_pair(matches)
                        f1.write(f"{variety_code}\t{os.path.basename(r1)}\t{os.path.basename(r2)}\n")
                f2.write(f"{variety_code}\t{variety}\n")