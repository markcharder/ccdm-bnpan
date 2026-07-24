#!/usr/bin/env python3

import argparse

def extract_fasta(
    fasta_file: str
    ):

    with open(fasta_file, 'r') as f:
        coords = dict()

        for line in f:

            if line.startswith('>'):
                header = line.strip()
                coords[header] = []

            else:
                coords[header] = \
                coords[header] + [line.strip().upper().split("")]
    
    return coords

def filter_coords(
    coords: dict,
    pad = 10
    ):

    coords_diffs = []

    for seq in coords.values():

