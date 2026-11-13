# Revisions
These are the steps taken to address reviewers' comments on the manuscript, which
was uploaded as a preprint at: https://doi.org/10.64898/2026.07.29.740455
Unfortunately, the variants used as input were not backed up, so the graph mapping
and variant calling pipeline was re-run. Two additional scripts are included for
downloading the reads from an aws bucket, where they are stored. Note, these reads
are due to be released from NCBI on 1st December under BioProject PRJNA1504441.
These are the scripts run on the HPC to download the reads:
```
scripts/slurm/dnldreads.slurm
scripts/slurm/dnldreads.sub
```
Simply run ```bash scripts/slurm/dnldreads.sub``` to download the reads to
```data/graphmap/reads_illumina```. Then, make the variety IDs config file for the
nextflow pipeline:
```bash
mkdir -p data/graphmap/conf/
python3 scripts/util/make_read_pairs.py
```
Then, run both graphmap and graphcall with:
```bash
sbatch scripts/slurm/graphvar.slurm
```

# Variants used as input for haplotype and selection analyses
These are extracted from `envs/bash/variants.sh` and `doc/command-line/VARIANTS.md`.
```bash

# Data directory variables:
WF_DIR="$WORKFLOWS/ccdm-bnpan/"
VARIANTS_DIR="$WF_DIR/data/variants"
MERGE_DIR="$VARIANTS_DIR/merged/"

# Variant file variables:
merge_full_bcf="$MERGE_DIR/merge_full.bcf.gz"
merge_filt_bcf="$MERGE_DIR/merge_filt.bcf.gz"
merge_filtsub_bcf="$MERGE_DIR/merge_filtsub.bcf.gz"
merge_basnpsmaf_bcf="$MERGE_DIR/merge_basnpsmaf.bcf.gz"

# Filter alleles to get biallelic, < 10% missing and >= 0.05 MAF and index.
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --include 'F_MISSING < 0.1 && MAF >= 0.05' \
        --min-alleles 2 --max-alleles 2 \
        $merge_full_bcf > $merge_filt_bcf

bcftools \
    index \
        --force \
        $merge_filt_bcf

# Get genotype info from the filtered BCF.
bcftools \
    +smpl-stats \
    $merge_filt_bcf > $MERGE_DIR/gt_info.txt

# Calculate number of variant sites overall.
nrecords=$(bcftools index --nrecords $merge_filt_bcf.csi)

# Get missing rate for each sample.
awk '$1 == "FLT0"' $MERGE_DIR/gt_info.txt | \
    cut -f 2,12 | awk -v nrec=$nrecords '{print $1"\t"$2/nrec}' \
    > $MERGE_DIR/gt_missrate.txt

# Keep good samples with <= 26% missing.
awk '$2 <= 0.26' $MERGE_DIR/gt_missrate.txt | \
    cut -f 1 > $MERGE_DIR/good_samples.txt

# Subset the filtered BCF to better quality samples.
bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --samples-file "$MERGE_DIR/good_samples.txt" \
        --trim-alt-alleles \
        $merge_filt_bcf | \
bcftools \
    reheader \
        --threads 30 \
        --samples-file $MERGE_DIR/good_samples.txt | \
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --min-alleles 2 --max-alleles 2 \
        # Keep only variants with <= 5% missing and >= 5% MAF after subset.
        --include 'F_MISSING <= 0.05 && MAF >= 0.05' \
        > $merge_filtsub_bcf

# Get only biallelic SNPs from the biallelic variants bcf.
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --include 'TYPE="snp" && strlen(REF)==1 && strlen(ALT)==1' \
        $merge_filtsub_bcf > $merge_basnpsmaf_bcf

```

# Fixing the allele frequency selection scan
The problem with the previous test was that the number of valid pairs changes depending on missing data (although all genotypes have a minimum of 95% call rate). This is a revised version of the function to perform the selection scan with more detailed comments.
```bash

normalise_deltap() {

    base_freqs="$1" # A TSV file containing variant positions and base allele fequencies.
    varpos="$2" # A TSV file containing variant positions and genotypes ('%CHROM\t%POS[\t%GT]\n').

    # Pass the base frequency file to awk as a variable.
    awk -v base_freqs="$base_freqs" '

    # Do this before parsing the varpos file.
    BEGIN {

        # Set field separator for input and output.
        FS = "\t";
        OFS = "\t";

        # This is hardcoded - see calculation above.
        sum_G_offdiag = 317.8384 # Set to sum of off-diagonals of kinship matrix.
        n_total = 20 # Number of individuals in the kinship matrix.
        total_pairs = n_total * (n_total - 1) # The total non-self pairs (number offdiags).
        theta_bar = sum_G_offdiag / total_pairs # The average relatedness for a pair.

        # While getline invoked on the base freqs file returns success (1).
        while ((getline < base_freqs) > 0) {

            # Add the base AF ($3) to an associative array, $1 = chrom, $2 = pos.
            p_map[$1"_"$2] = $3;

        }

        close(base_freqs);

        # Print the header for the output.
        print "CHROM", "POS", "LOCAL_SHARING", "EXPECTED_SHARING", "DEVIATION", "VALID_N", "P_BASE", "P_MODERN";

    }

    # This applies to varpos.
    {
        # Chromosome, position and p - the base allele frequency.
        chr = $1;
        pos = $2;
        p = p_map[chr"_"pos];
    
        # Skip of base AF is missing or 0.
        if (p == "" || p == 0) next; 

        # Set counter variables to 0.
        sum_c = 0;
        sum_c_sq = 0;
        valid_n = 0;
        sum_raw_alleles = 0

        # From 3 to NF captures genotypes.
        for (i = 3; i <= NF; i++) {

            # Skip the genotype if it is missing.
            if ($i ~ /\.\/\./ || $i ~ /\./) continue;

            # Get the alleles in the genotype.
            split($i, alleles, "/|\\|");

            # The genotype score is either 0 or 2 (0 + 0 or 1 + 1).
            genotype_score = alleles[1] + alleles[2];

            # Add to this to calculate p_modern.
            sum_raw_alleles += genotype_score;

            # Subtract 2 x base AF (p) to make scale the genotype to the mean base
            # AF.
            c = genotype_score - (2 * p);
            sum_c += c;
            sum_c_sq += (c * c);

            # Increment valid number of samples seen.
            valid_n++;

        }

        # Skip if we have < 2 valid samples (should not happen).
        if (valid_n < 2) next;

        # Non-self square (offdiag) comparisons.
        pairs = valid_n * (valid_n - 1);

        # Local sharing is the sum of centred GTs - sum of c^2.
        local_sharing = (sum_c * sum_c) - sum_c_sq;

        # Expected sharing NOW SCALES AVERAGE KINSHIP WITH VALID PAIRS.
        expected_sharing = pairs * theta_bar * (2 * p * (1 - p));

        # Deviation is (local - expected ) / pairs. The average AF deviation per pair.
        deviation = (local_sharing - expected_sharing) / pairs;

        p_modern =  sum_raw_alleles / (2 * valid_n)
        print chr, pos, local_sharing, expected_sharing, deviation, valid_n, p, p_modern;

    }' $varpos

}

export -f normalise_deltap

```
