#!/usr/bin/env python3

import argparse
import sys

parser = argparse.ArgumentParser(
    description=(
        "Split a RagTag-scaffolded FASTA into A-genome, C-genome, "
        "and unplaced files using a TSV that maps reference scaffold "
        "names to genome (A or C)."
    )
)

parser.add_argument(
    "--fasta",
    "-f",
    required=True,
    help="Input scaffolded FASTA file (ragtag.scaffold.fasta)."
)

parser.add_argument(
    "--tsv",
    "-t",
    required=True,
    help=(
        "Two-column TSV mapping reference scaffold names to genome. "
        "Column 1: scaffold name (e.g. ChrA01). "
        "Column 2: genome label (A or C). No header row."
    )
)

parser.add_argument(
    "--sample",
    "-s",
    required=True,
    help="Sample name used as prefix for output file names."
)

parser.add_argument(
    "--ragtag-suffix",
    default="_RagTag",
    help="Suffix appended by RagTag to placed scaffold names (default: _RagTag)."
)


def load_tsv(tsv_file: str) -> dict:
    mapping = {}
    with open(tsv_file, "r") as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                sys.stderr.write(
                    f"WARNING: skipping malformed line {lineno} in TSV: {line!r}\n"
                )
                continue
            scaffold, genome = parts
            genome = genome.upper()
            if genome not in ("A", "C"):
                sys.stderr.write(
                    f"WARNING: unrecognised genome label {genome!r} on line {lineno}, skipping.\n"
                )
                continue
            mapping[scaffold] = genome
    return mapping


def split_fasta(fasta_file: str, mapping: dict, sample: str, ragtag_suffix: str):
    out_paths = {
        "A":        f"{sample}_A.fasta",
        "C":        f"{sample}_C.fasta",
        "unplaced": f"{sample}_unplaced.fasta",
    }

    handles = {k: open(v, "w") for k, v in out_paths.items()}

    try:
        current_handle = None

        with open(fasta_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    # Extract bare sequence name (no ">" or description)
                    seq_name = line[1:].split()[0]

                    # Strip RagTag suffix to recover the reference scaffold name
                    if seq_name.endswith(ragtag_suffix):
                        ref_name = seq_name[: -len(ragtag_suffix)]
                        genome = mapping.get(ref_name)
                    else:
                        genome = None

                    current_handle = handles[genome if genome else "unplaced"]
                    current_handle.write(line)
                else:
                    if current_handle is None:
                        sys.stderr.write("Malformed FASTA: sequence line before any header.\n")
                        sys.exit(1)
                    current_handle.write(line)
    finally:
        for h in handles.values():
            h.close()

    for label, path in out_paths.items():
        sys.stderr.write(f"Written {label}: {path}\n")


if __name__ == "__main__":
    args = parser.parse_args()
    mapping = load_tsv(args.tsv)
    split_fasta(
        fasta_file=args.fasta,
        mapping=mapping,
        sample=args.sample,
        ragtag_suffix=args.ragtag_suffix,
    )
