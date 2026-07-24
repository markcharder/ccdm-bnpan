#!/usr/bin/env python3

import sys
import numpy as np

mash_file = sys.argv[1]

def analyse_mashdist(filepath):
    same, diff = [], []
    names = []
    
    with open(filepath) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue  # Skip header line
            parts = line.strip().split('\t')
            genome1 = parts[0].split('#')[2]
            names.append(parts[0])
            
            for j, dist in enumerate(parts[1:], 1):
                if dist and float(dist) > 0:
                    genome2 = names[i-j].split('#')[2]
                    (same if genome1 == genome2 else diff).append(float(dist))
    
    for label, data in [("Same genome", same), ("Different genomes", diff)]:
        print(f"{label}:")
        print(f"  Min: {np.min(data):.6f}  Max: {np.max(data):.6f}")
        print(f"  Mean: {np.mean(data):.6f}  Median: {np.median(data):.6f}\n")

analyse_mashdist(mash_file)