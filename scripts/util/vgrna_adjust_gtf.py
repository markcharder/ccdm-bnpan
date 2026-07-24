#!/usr/bin/env python3

import argparse
import re
import subprocess
import sys

from bisect import bisect_right
from collections import defaultdict

def get_contig_offsets(offsets_path: str, gbz: bool) -> dict:

    if gbz:
        # Run vg paths to get all paths and lengths.
        sys.stderr.write(f"Offsets in gbz. Running vg paths...\n")
        result = subprocess.run(
            ["vg", "paths", "-x", offsets_path, "--lengths"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True,
        )
        offsets_file = result.stdout.splitlines()
        sys.stderr.write(f"Got {len(offsets_file)} paths from vg.\n")
    else:
        sys.stderr.write(f"Offsets in text file. Reading from {offsets_path}...\n")
        with open(offsets_path, "r") as f:
            offsets_file = f.read().splitlines()
        sys.stderr.write(f"Got {len(offsets_file)} paths from text file.\n")

    # Nested default dictionary to store offsets.
    offsets = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    sys.stderr.write(f"Parsing offsets...\n")
    # For each line in the output of vg.
    for line in offsets_file:

        try:
            # Get the path name and length.
            path_name, path_length = line.strip().split()

            # The sample, haplotype and contig are in the pan-SN spec name.
            parts = path_name.split("#")
            sample, haplotype, contig = parts[0], parts[1], parts[2]

            # The offset is the 4th part if this is not a reference path.
            # Otherwise it is set to zero (if it is reference).
            offset = int(parts[3]) if len(parts) == 4 else 0

            # Provide offsets and path lengths for finding intervals in contigs.
            offsets[sample][haplotype][contig].append((offset, int(path_length)))

        except (IndexError, ValueError):
            sys.stderr.write(f"WARNING: malformed line: {line}\n")

    sys.stderr.write(f"Parsed offsets for {len(offsets)} samples.\n")
    sys.stderr.write(f"Sorting offsets...\n")

    # Sort each contig's fragment list by offset ascending.
    for sample in offsets:

        for haplotype in offsets[sample]:

            for contig in offsets[sample][haplotype]:
                offsets[sample][haplotype][contig].sort()
    sys.stderr.write(f"Offsets sorted.\n")

    return offsets

def find_fragment(fragments: list, feat_start: int, feat_end: int):
    # Starts is a list of fragment offsets.
    starts = [f[0] for f in fragments]

    # bisect_right returns the index where feat_start would be
    # inserted to maintain sorted order.
    i = bisect_right(starts, feat_start) - 1

    # If the feature start is before the first fragment, we can't place it.
    if i < 0:
        return None, None, f"Feature start {feat_start} is before all fragments"

    frag_offset, frag_length = fragments[i]
    frag_end = frag_offset + frag_length

    warning = None

    # If the feature extends beyond the end of the fragment, we can't place it.
    if feat_end > frag_end:
        
        if feat_start >= frag_end:
            return None, None, f"Feature [{feat_start}, {feat_end}] is after fragment end at {frag_end}"

        else:
            return None, None, f"Feature [{feat_start}, {feat_end}] spans fragment boundary at {frag_end}"

    return frag_offset, frag_offset, None

def adjust_gtf(
    gtf_path: str,
    offsets: dict,
    sample: str,
    haplotype: str,
    verbose: bool = False
):
    # Get the offsets for this sample and haplotype.
    hap_offsets = offsets.get(sample, {}).get(haplotype, {})
    pansn_prefix = f"{sample}#{haplotype}#"

    # Open the GTF and read by line.
    with open(gtf_path, "r") as f:

        # For each line, adjust coordinates, if possible.
        for line in f:

            # Print comment lines.
            if line.startswith("#"):
                print(line.strip())
                continue

            # Split to get the feature coordinates etc.
            parts = line.strip().split("\t")

            # Ignore malformed lines.
            if len(parts) < 5:
                print(line.strip())
                continue

            # Get the contig name.
            contig = parts[0]

            # Strip PanSN prefix if present to get bare contig name.
            if contig.startswith(pansn_prefix):
                contig_id = contig.split("#")[2]

            else:
                contig_id = contig

            # GTF is 1-based, subtract 1 to get 0-based coordinates for search.
            feat_start = int(parts[3]) - 1
            feat_end   = int(parts[4])      # End is inclusive in GTF.

            # Get the fragments for this contig.
            fragments = hap_offsets.get(contig_id)

            # If there are no fragments for this contig, write a warning and print the line.
            if not fragments:

                if verbose:
                    sys.stderr.write(f"WARNING: no fragments found for contig '{contig_id}'\n")

                print(line.strip())
                continue

            # Find the fragment offset and whether the feature fits within
            # a path.
            frag_offset, phase_block, warning = find_fragment(fragments, feat_start, feat_end)

            if warning:

                if verbose:
                    sys.stderr.write(f"WARNING: {contig_id}: {warning}\n")

                print(line.strip())
                continue

            # Rename contig to full PanSN path name and convert back to 1-based.
            parts[0] = f"{sample}#{haplotype}#{contig_id}#{phase_block}"
            parts[3] = str(feat_start - frag_offset + 1)
            parts[4] = str(feat_end   - frag_offset)

            print("\t".join(parts))

def main():
    parser = argparse.ArgumentParser(
        description = "Adjust GTF coordinates to match offsets in a minigraph cactus GBZ."
    )
    parser.add_argument(
        "--gbz",
        default = False,
        help = "Whether or not the offsets are in a gbz file.",
        action = "store_true"
    )
    parser.add_argument(
        "--verbose",
        default = False,
        help = "Whether or not to print warnings about features that can't be placed.",
        action = "store_true"
    )
    parser.add_argument(
        "offsets",
        help = "Path to the offsets file. Can either be a gbz or output of vg paths --lengths"
    )
    parser.add_argument(
        "gtf",
        help = "Path to the GTF file to adjust."
    )
    parser.add_argument(
        "sample",
        help = "Sample name."
    )
    parser.add_argument(
        "haplotype",
        help = "Haplotype name."
    )
    args = parser.parse_args()

    offsets = get_contig_offsets(args.offsets, args.gbz)
    adjust_gtf(args.gtf, offsets, args.sample, args.haplotype, args.verbose)

if __name__ == "__main__":
    main()
