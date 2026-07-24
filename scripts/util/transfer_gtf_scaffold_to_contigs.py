#!/usr/bin/env python3

import argparse
import sys
from bisect import bisect_right

def parse_paf(
    paf_file: str
    ):

    scaffold_to_contig = dict()
    track_primaries = dict()

    with open(paf_file, 'r') as f:

        for line in f:
            parts = line.strip().split('\t')
            contig_name = parts[0]
            contig_length = int(parts[1])
            contig_start = int(parts[2])
            contig_end = int(parts[3])
            contig_strand = parts[4]
            scaffold_name = parts[5]
            scaffold_start = int(parts[7])
            scaffold_end = int(parts[8])
            residue_matches = int(parts[9])
            block_len = int(parts[10])

            # Parse optional tags by name rather than fixed index
            tags = {p.split(':')[0]: p for p in parts[12:]}
            is_primary = tags.get('tp') == 'tp:A:P'

            # Exact match: every aligned base is a match
            is_exact = block_len > 0 and residue_matches == block_len

            if is_exact and is_primary and block_len == contig_length:

                if scaffold_name not in scaffold_to_contig:
                    scaffold_to_contig[scaffold_name] = []

                if contig_name not in track_primaries:
                    track_primaries[contig_name] = 0

                scaffold_to_contig[scaffold_name].append(
                    (
                        scaffold_start, scaffold_end,
                        contig_start, contig_end,
                        contig_strand, contig_name
                    )
                )
                track_primaries[contig_name] += 1

    # Sort each scaffold's alignments by scaffold_start once
    for alns in scaffold_to_contig.values():
        alns.sort(key=lambda a: a[0])

    # Precompute starts lists once per scaffold
    scaffold_starts = {
        scaf: [a[0] for a in alns]
        for scaf, alns in scaffold_to_contig.items()
    }

    return scaffold_to_contig, track_primaries, scaffold_starts

def get_skip_transcripts(
    gtf_file: str,
    scaffold_to_contig: dict,
    scaffold_starts: dict
    ):
    """
    First pass: scan transcript lines and collect transcript_ids whose full
    span is not contained within a single contig. These will be skipped in
    the second pass.
    """
    skip_transcripts = set()

    def extract_transcript_id(attrs):
        for attr in attrs.split(';'):
            attr = attr.strip()
            if attr.startswith('transcript_id'):
                return attr.split('"')[1] if '"' in attr else attr.split()[-1]
        return None

    with open(gtf_file, 'r') as f:

        for line in f:

            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')

            if len(parts) < 9 or parts[2] != 'transcript':
                continue

            scaffold_name = parts[0]
            feature_start = int(parts[3])
            feature_end = int(parts[4])
            transcript_id = extract_transcript_id(parts[8]) or parts[8].strip()

            if transcript_id is None:
                continue

            if scaffold_name not in scaffold_to_contig:
                print(
                    f"WARNING: transcript '{transcript_id}' skipped — scaffold "
                    f"'{scaffold_name}' not in PAF.",
                    file=sys.stderr,
                )
                skip_transcripts.add(transcript_id)
                continue

            upper = bisect_right(scaffold_starts[scaffold_name], feature_start)
            found = any(
                feature_start >= s and feature_end <= e
                for s, e, *_ in scaffold_to_contig[scaffold_name][:upper]
            )

            if not found:
                print(
                    f"WARNING: transcript '{transcript_id}' skipped — span "
                    f"{scaffold_name}:{feature_start}-{feature_end} not contained in a single contig.",
                    file=sys.stderr,
                )
                skip_transcripts.add(transcript_id)

    return skip_transcripts


def parse_gtf(
    gtf_file: str,
    scaffold_to_contig: dict,
    track_primaries: dict,
    scaffold_starts: dict,
    skip_transcripts: set
    ):

    def extract_transcript_id(attrs):
        for attr in attrs.split(';'):
            attr = attr.strip()
            if attr.startswith('transcript_id'):
                return attr.split('"')[1] if '"' in attr else attr.split()[-1]
        return None

    with open(gtf_file, 'r') as f:

        for line in f:

            if line == "\n" or line.startswith('#'):
                print(line.strip())
                continue

            else:
                parts = line.strip().split('\t')
                transcript_id = extract_transcript_id(parts[8]) if len(parts) >= 9 else None

                if transcript_id in skip_transcripts:
                    continue

                scaffold_name = parts[0]
                feature_start = int(parts[3])
                feature_end = int(parts[4])
                feature_strand = parts[6]

                if scaffold_name not in scaffold_to_contig:
                    print(
                        f"WARNING: scaffold '{scaffold_name}' not in PAF — skipping feature "
                        f"{feature_start}-{feature_end}.",
                        file=sys.stderr,
                    )
                    continue

                upper = bisect_right(scaffold_starts[scaffold_name], feature_start)

                best = None
                palign_len = 0

                for scaffold_start, scaffold_end, contig_start, contig_end, contig_strand, contig_name \
                    in scaffold_to_contig[scaffold_name][:upper]:

                    if feature_start >= scaffold_start and feature_end <= scaffold_end:
                        align_len = scaffold_end - scaffold_start

                        if track_primaries[contig_name] > 1:
                            print(
                                f"WARNING: contig '{contig_name}' has multiple primary alignments "
                                f"— skipping feature {scaffold_name}:{feature_start}-{feature_end}.",
                                file=sys.stderr
                            )
                            continue

                        if align_len > palign_len:

                            if contig_strand == '+':
                                new_start = feature_start - scaffold_start
                                new_end = feature_end - scaffold_start
                                new_strand = feature_strand

                            else:
                                new_start = contig_end - (feature_end - scaffold_start) + 1
                                new_end = contig_end - (feature_start - scaffold_start) + 1

                                if feature_strand == '+':
                                    new_strand = '-'

                                elif feature_strand == '-':
                                    new_strand = '+'

                                else:
                                    new_strand = feature_strand

                            best = (contig_name, new_start, new_end, new_strand)
                            palign_len = align_len

                if best is not None:
                    parts[0], parts[3], parts[4], parts[6] = best[0], str(best[1]), str(best[2]), best[3]
                    print('\t'.join(parts))

def read_fasta(fasta_file: str) -> dict:
    """
    Returns dict: sequence_name -> sequence string (upper case, no newlines).
    """
    seqs = {}
    name = None
    buf = []

    with open(fasta_file, 'r') as f:

        for line in f:
            line = line.strip()

            if line.startswith('>'):

                if name is not None:
                    seqs[name] = ''.join(buf).upper()
                name = line[1:].split()[0]
                buf = []

            else:
                buf.append(line)

    if name is not None:
        seqs[name] = ''.join(buf).upper()

    return seqs

REV_COMP = str.maketrans('ACGTacgt', 'TGCAtgca')

def revcomp(seq: str) -> str:
    return seq.translate(REV_COMP)[::-1]


def extract_cds_from_gtf(output_gtf: str, contig_seqs: dict) -> dict:
    """
    Parses the transferred GTF and extracts the concatenated CDS sequence
    for each transcript, using the contig sequences.

    Returns dict: transcript_id -> extracted sequence string.
    """
    from collections import defaultdict

    # transcript_id -> list of (start_1based, end_1based, strand, contig_name)
    cds_features = defaultdict(list)

    with open(output_gtf, 'r') as f:

        for line in f:

            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')

            if len(parts) < 9 or parts[2] != 'CDS':
                continue

            contig_name = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attrs = parts[8]

            # Extract transcript_id from the attributes field
            transcript_id = None

            for attr in attrs.split(';'):
                attr = attr.strip()

                if attr.startswith('transcript_id'):
                    transcript_id = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                    break

            if transcript_id is None:
                continue

            cds_features[transcript_id].append((start, end, strand, contig_name))

    extracted = {}

    for transcript_id, exons in cds_features.items():
        strand = exons[0][2]

        # Sort by coordinate; reverse order for minus-strand transcripts
        exons_sorted = sorted(exons, key=lambda e: e[0], reverse=(strand == '-'))
        seq_parts = []

        for start, end, strand, contig_name in exons_sorted:

            if contig_name not in contig_seqs:
                print(
                    f"WARNING: contig '{contig_name}' not found in FASTA — "
                    f"skipping transcript '{transcript_id}'.",
                    file=sys.stderr,
                )
                seq_parts = []
                break

            # GTF is 1-based closed; Python slicing is 0-based half-open
            exon_seq = contig_seqs[contig_name][start - 1:end]

            if strand == '-':
                exon_seq = revcomp(exon_seq)
            seq_parts.append(exon_seq)

        if seq_parts:
            extracted[transcript_id] = ''.join(seq_parts)

    return extracted

def test_transfer(output_gtf: str, contigs_file: str, cds_file: str):
    """
    Compares CDS sequences extracted from the transferred GTF + contig FASTA
    against the BRAKER3 coding sequence FASTA.
    """
    print("Running transfer test...", file=sys.stderr)

    contig_seqs = read_fasta(contigs_file)
    braker_cds  = read_fasta(cds_file)
    extracted   = extract_cds_from_gtf(output_gtf, contig_seqs)

    n_match   = 0
    n_mismatch = 0
    n_missing  = 0

    for transcript_id, ref_seq in braker_cds.items():

        if transcript_id not in extracted:
            print(
                f"MISSING: '{transcript_id}' not found in transferred GTF.",
                file=sys.stderr,
            )
            n_missing += 1
            continue

        if extracted[transcript_id] == ref_seq:
            n_match += 1

        else:

            if ref_seq[-1] == "N":
                
                if extracted[transcript_id][:-1] == ref_seq[:-1]:
                    print(
                        f"MATCH with trailing N: '{transcript_id}' — "
                        f"extracted length {len(extracted[transcript_id])}, "
                        f"reference length {len(ref_seq)}.",
                        f"extracted start: {extracted[transcript_id][0:10]}",
                        f"reference start: {ref_seq[0:10]}",
                        f"extracted end: {extracted[transcript_id][-10:]}",
                        f"reference end: {ref_seq[-10:]}",
                        file=sys.stderr,
                    )
                    n_match += 1
                    continue

            elif ref_seq[0] == "N":

                if extracted[transcript_id][1:] == ref_seq[1:]:
                    print(
                        f"MATCH with leading N: '{transcript_id}' — "
                        f"extracted length {len(extracted[transcript_id])}, "
                        f"reference length {len(ref_seq)}.",
                        f"extracted start: {extracted[transcript_id][0:10]}",
                        f"reference start: {ref_seq[0:10]}",
                        f"extracted end: {extracted[transcript_id][-10:]}",
                        f"reference end: {ref_seq[-10:]}",
                        file=sys.stderr,
                    )
                    n_match += 1
                    continue

            else:
                print(
                    f"MISMATCH: '{transcript_id}' — "
                    f"extracted length {len(extracted[transcript_id])}, "
                    f"reference length {len(ref_seq)}.",
                    f"extracted start: {extracted[transcript_id][0:10]}",
                    f"reference start: {ref_seq[0:10]}",
                    f"extracted end: {extracted[transcript_id][-10:]}",
                    f"reference end: {ref_seq[-10:]}",
                    file=sys.stderr,
                )
                n_mismatch += 1

    print(
        f"Test complete.  Match: {n_match}  Mismatch: {n_mismatch}  Missing: {n_missing}",
        file=sys.stderr,
    )

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Transfer GTF annotations from scaffold coordinates to contig "
            "coordinates using a PAF alignment file (contigs → scaffolds)."
        )
    )
    p.add_argument("--paf", required=True,
                   help="PAF file: contigs (query) aligned to scaffolds (target)")
    p.add_argument("--gtf", required=True,
                   help="Input GTF with annotations on scaffold sequences")
    p.add_argument("--output", required=True,
                   help="Output GTF with annotations transferred to contig sequences")
    p.add_argument("--test", action="store_true",
                   help="Run a test comparing extracted CDS to BRAKER3 sequences")
    p.add_argument("--contigs",
                   help="Contig FASTA file (required if --test is set)")
    p.add_argument("--cds",
                   help="BRAKER3 coding sequence FASTA file (required if --test is set)")
    return p.parse_args()


def main():
    args = parse_args()

    if args.test and (not args.contigs or not args.cds):
        print("ERROR: --test requires --contigs and --cds.", file=sys.stderr)
        sys.exit(1)

    scaffold_to_contig, track_primaries, scaffold_starts = parse_paf(args.paf)
    skip_transcripts = get_skip_transcripts(args.gtf, scaffold_to_contig, scaffold_starts)

    original_stdout = sys.stdout
    with open(args.output, 'w') as fh:
        sys.stdout = fh
        parse_gtf(args.gtf, scaffold_to_contig, track_primaries, scaffold_starts, skip_transcripts)
    sys.stdout = original_stdout

    if args.test:
        test_transfer(args.output, args.contigs, args.cds)


if __name__ == '__main__':
    main()