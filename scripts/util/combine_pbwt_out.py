#!/usr/bin/env python3

import argparse

# Read in samples for reference and target individuals
# and assign numerical 0-based genotypes to them.
def get_genotypes(ref_samples, target_samples):
    ref_genotypes = dict()
    target_genotypes = dict()

    with open(ref_samples, 'r') as f:
        idx = 0
        for line in f:
            sample_id = line.strip()
            if sample_id:
                ref_genotypes[sample_id] = idx
                idx += 2

    with open(target_samples, 'r') as f:
        idx = 0
        for line in f:
            sample_id = line.strip()
            if sample_id:
                target_genotypes[sample_id] = idx
                idx += 2

    return ref_genotypes, target_genotypes

# Read SNP positions and chromosomes from sites file.
# Returns a list where index = SNP index (0-based).
# Each element is a tuple: (chromosome, base_pair_position)
def get_snp_positions(sites_file):
    snp_positions = []
    
    with open(sites_file, 'r') as f:

        for line in f:
            fields = line.strip().split('\t')

            if len(fields) >= 2:
                chromosome = fields[0]
                bp_position = int(fields[1])
                snp_positions.append((chromosome, bp_position))

    return snp_positions

# Read matches file and filter to one genotype per sample.
# Returns a list of matches, where each match is a dict with:
#   target_geno, ref_geno, start_snp_idx, end_snp_idx, num_snps
def read_matches(matches_file):
    matches = []
    
    with open(matches_file, 'r') as f:

        for line in f:
            fields = line.strip().split('\t')

            if len(fields) >= 6 and fields[0] == 'MATCH':
                target_geno = int(fields[1])
                ref_geno = int(fields[2])
                start_snp_idx = int(fields[3])
                end_snp_idx = int(fields[4]) - 1
                num_snps = int(fields[5])

                # Only keep even-numbered genotypes (one per sample)
                # This filters out the redundant genotype for each sample
                if target_geno % 2 == 0 and ref_geno % 2 == 0:
                    matches.append({
                        'target_geno': target_geno,
                        'ref_geno': ref_geno,
                        'start_snp_idx': start_snp_idx,
                        'end_snp_idx': end_snp_idx,
                        'num_snps': num_snps
                    })
    
    return matches

# Convert SNP-based matches to base pair coordinates.
# Takes matches list and snp_positions list.
# Returns matches with added fields: chromosome, start_bp, end_bp, bp_length
def convert_matches_to_bp(matches, snp_positions):
    converted_matches = []

    for match in matches:
        start_idx = match['start_snp_idx']
        end_idx = match['end_snp_idx']

        # Get chromosome and base pair positions for start and end SNPs
        chrom_start, bp_start = snp_positions[start_idx]
        chrom_end, bp_end = snp_positions[end_idx]

        # Sanity check: start and end should be on same chromosome
        if chrom_start != chrom_end:
            continue  # Skip matches spanning multiple chromosomes

        # Add base pair coordinate information to the match
        match['chromosome'] = chrom_start
        match['start_bp'] = bp_start
        match['end_bp'] = bp_end
        match['bp_length'] = bp_end - bp_start + 1  # Inclusive range

        converted_matches.append(match)

    return converted_matches

# Identify set maximal matches.
# A match is set maximal if it's the ONLY longest match ending at that SNP position
# FOR THAT TARGET SAMPLE. Each target sample has its own set maximals.
# Takes matches list and adds a 'is_set_maximal' boolean field to each match.
# Returns the modified matches list.
def identify_set_maximals(matches):

    # Group matches by (target_geno, end_snp_idx) - set maximals are per target sample.
    matches_by_target_and_end = {}

    for match in matches:
        target_geno = match['target_geno']
        end_idx = match['end_snp_idx']
        key = (target_geno, end_idx)

        if key not in matches_by_target_and_end:
            matches_by_target_and_end[key] = []

        matches_by_target_and_end[key].append(match)
    
    # For each (target, end) combination, find the longest match(es)
    for key, end_matches in matches_by_target_and_end.items():

        # Find maximum length among all matches for this target ending at this position
        max_length = max(m['num_snps'] for m in end_matches)
        
        # Find all matches with that maximum length
        longest_matches = [m for m in end_matches if m['num_snps'] == max_length]
        
        # Mark as set maximal only if there's exactly ONE longest match for this target
        for match in end_matches:

            if len(longest_matches) == 1 and match in longest_matches:
                match['is_set_maximal'] = True

            else:
                match['is_set_maximal'] = False

    return matches

# Analyze match coverage efficiently using an event-based approach.
# Outputs a TSV file with coverage information for each contiguous region per target sample.
# Also outputs a second file with unique set maximal coordinates for single set maximal regions.
def analyze_and_output_coverage(matches, ref_genotypes, target_genotypes, output_file, set_maximal_file):
    # Reverse the genotype dictionaries to map genotype ID -> sample name
    ref_geno_to_sample = {v: k for k, v in ref_genotypes.items()}
    target_geno_to_sample = {v: k for k, v in target_genotypes.items()}
    
    # Group matches by target sample, then by chromosome
    matches_by_target_and_chrom = {}

    for match in matches:
        target_geno = match['target_geno']
        chrom = match['chromosome']

        if target_geno not in matches_by_target_and_chrom:
            matches_by_target_and_chrom[target_geno] = {}
        
        if chrom not in matches_by_target_and_chrom[target_geno]:
            matches_by_target_and_chrom[target_geno][chrom] = []

        matches_by_target_and_chrom[target_geno][chrom].append(match)
    
    # Open output files
    with open(output_file, 'w') as out, open(set_maximal_file, 'w') as sm_out:

        # Write headers
        out.write('\t'.join([
            'target_sample',
            'chromosome',
            'start_bp',
            'end_bp',
            'n_maximal_matches',
            'n_set_maximal',
            'maximal_ids',
            'set_maximal_ids'
        ]) + '\n')
        
        sm_out.write('\t'.join([
            'target_sample',
            'chromosome',
            'region_start_bp',
            'region_end_bp',
            'ref_sample',
            'match_start_bp',
            'match_end_bp'
        ]) + '\n')

        # Process each target sample
        for target_geno in sorted(matches_by_target_and_chrom.keys()):
            target_sample = target_geno_to_sample.get(target_geno, f"geno_{target_geno}")
            
            # Process each chromosome for this target sample
            for chrom in sorted(matches_by_target_and_chrom[target_geno].keys()):
                chrom_matches = matches_by_target_and_chrom[target_geno][chrom]
                
                # Create events: (position, event_type, match_object)
                # event_type: 1 for start, -1 for end
                events = []

                for match in chrom_matches:
                    events.append((match['start_bp'], 1, match))
                    events.append((match['end_bp'] + 1, -1, match))  # +1 because end is inclusive
                
                # Sort events by position, with starts before ends at same position
                events.sort(key=lambda x: (x[0], -x[1]))

                # Track active matches
                active_matches = []
                active_set_maximals = []
                prev_pos = None
                
                for pos, event_type, match in events:

                    # Output region for previous position if coverage changed
                    if prev_pos is not None and pos > prev_pos and len(active_matches) > 0:

                        # Format match information
                        maximal_info = []

                        for m in active_matches:
                            ref_sample = ref_geno_to_sample.get(m['ref_geno'], f"geno_{m['ref_geno']}")
                            maximal_info.append(f"{ref_sample}:{m['start_bp']}-{m['end_bp']}")
                        
                        set_maximal_info = []
                        for m in active_set_maximals:
                            ref_sample = \
                                ref_geno_to_sample.get(
                                    m['ref_geno'], f"geno_{m['ref_geno']}"
                                )
                            set_maximal_info.append(
                                f"{ref_sample}:{m['start_bp']}-{m['end_bp']}"
                            )
                        
                        # Write output to main file
                        out.write('\t'.join([
                            target_sample,
                            chrom,
                            str(prev_pos),
                            str(pos - 1),
                            str(len(active_matches)),
                            str(len(active_set_maximals)),
                            ';'.join(maximal_info) if maximal_info else '.',
                            ';'.join(set_maximal_info) if set_maximal_info else '.'
                        ]) + '\n')
                        
                        # Write to set maximal file if exactly one set maximal
                        if len(active_set_maximals) == 1:
                            sm = active_set_maximals[0]
                            ref_sample = \
                                ref_geno_to_sample.get(
                                    sm['ref_geno'], f"geno_{sm['ref_geno']}"
                                    )
                            sm_out.write('\t'.join([
                                target_sample,
                                chrom,
                                str(prev_pos),
                                str(pos - 1),
                                ref_sample,
                                str(sm['start_bp']),
                                str(sm['end_bp'])
                            ]) + '\n')
                    
                    # Update active matches
                    if event_type == 1:  # Start event
                        active_matches.append(match)

                        if match['is_set_maximal']:
                            active_set_maximals.append(match)

                    else:  # End event
                        active_matches.remove(match)

                        if match['is_set_maximal']:
                            active_set_maximals.remove(match)
                    
                    prev_pos = pos

def main():
    parser = argparse.ArgumentParser(
        description = \
            "Identify high confidence ancestral haplotype blocks \
            based on pbwt maximal matches."
    )

    parser.add_argument(
        "--ref-samples", "-r",
        help = "List of reference samples, one per line. Order must match genotypes.",
        required = True
    )
    
    parser.add_argument(
        "--target-samples", "-t",
        help = "List of target samples, one per line. Order must match genotypes.",
        required = True
    )

    parser.add_argument(
        "--matches", "-m",
        help = "Maximal matches from pbwt. A TSV file with columns:\n\
            \t1) The string 'MATCH' (will be ignored).\n\
            \t2) Target genotype, 0-indexed numerical, 2 per sample.\n\
            \t3) Reference genotype, 0-indexed numerical, 2 per sample.\n\
            \t4) Match start SNP index 0-indexed.\n\
            \t5) Match end SNP index 0-indexed.\n\
            \t6) Number of SNPs in match.\n",
        required = True
    )

    parser.add_argument(
        "--sites", "-s",
        help = "The base pair positions of SNPs. A TSV file with columns:\n\
        \t1) Chromosome.\n\
        \t2) Base pair position of SNP (1-based).\n\
        \t3) SNP reference allele.\n\
        \t4) SNP alternate allele.\n",
        required = True
    )

    parser.add_argument(
        "--output", "-o",
        help = "Output TSV file with match coverage information.",
        required = True
    )

    parser.add_argument(
        "--set-maximal-output", "-sm",
        help = "Output TSV file with unique set maximal coordinates for single set maximal regions.",
        required = True
    )

    args = parser.parse_args()

    # Step 1: Read sample genotypes
    print("Reading sample genotypes...")
    ref_genotypes, target_genotypes = get_genotypes(args.ref_samples, args.target_samples)
    print(f"  Reference samples: {len(ref_genotypes)}")
    print(f"  Target samples: {len(target_genotypes)}")

    # Step 2: Read SNP positions
    print("Reading SNP positions...")
    snp_positions = get_snp_positions(args.sites)
    print(f"  Total SNPs: {len(snp_positions)}")

    # Step 3: Read and filter matches
    print("Reading matches...")
    matches = read_matches(args.matches)
    print(f"  Matches loaded: {len(matches)}")

    # Step 4: Convert matches to base pair coordinates
    print("Converting matches to base pair coordinates...")
    matches = convert_matches_to_bp(matches, snp_positions)
    print(f"  Matches with valid coordinates: {len(matches)}")

    # Step 5: Identify set maximal matches
    print("Identifying set maximal matches...")
    matches = identify_set_maximals(matches)
    set_maximal_count = sum(1 for m in matches if m['is_set_maximal'])
    print(f"  Set maximal matches: {set_maximal_count}")

    # Step 6: Analyze and output coverage
    print(f"Analyzing coverage and writing to {args.output}...")
    print(f"  Set maximal regions will be written to {args.set_maximal_output}")
    analyze_and_output_coverage(matches, ref_genotypes, target_genotypes, args.output, args.set_maximal_output)
    print("Done!")

if __name__ == "__main__":
    main()