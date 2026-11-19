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

# REVISION 1: kinship and selection analyses
For context surrounding revisions to the manuscript, we bring together code and intermediate files used for data preparation and subsequent processing and analysis. Comments marked with "# REVISION" are revisions made to code that impact results presented in the manuscript. These are applied in response to reviewer comments.

## Variants used for these analyses
Here, we assume earlier data preparation steps have been run so that the directory `data/variants/merged` has been created and it contains the full merged set of variants `merge_full.bcf.gz`. This is filtered as before and used to estimate kinship with no changes to the code.
```bash

# Data directory variables.
WF_DIR="$WORKFLOWS/ccdm-bnpan/"
VARIANTS_DIR="$WF_DIR/data/variants/"
MERGE_DIR="$VARIANTS_DIR/merged/"
SUBSET_DIR="$VARIANTS_DIR/subset/"

mkdir -p "$SUBSET_DIR"

# Variant file variables.
merge_full_bcf="$MERGE_DIR/merge_full.bcf.gz"
merge_filt_bcf="$MERGE_DIR/merge_filt.bcf.gz"
merge_filtsub_bcf="$MERGE_DIR/merge_filtsub.bcf.gz"
merge_ba_bcf="$MERGE_DIR/merge_ba.bcf.gz"
merge_basnpsmaf_bcf="$MERGE_DIR/merge_basnpsmaf.bcf.gz"
merge_basnpsmaf_vcf="$SUBSET_DIR/merge_basnpsmaf.vcf.gz"
merge_basnpsmaf_impute_vcf="$SUBSET_DIR/merge_basnpsmaf_impute"
merge_basnpsmaf_impute_bcf="$SUBSET_DIR/merge_basnpsmaf_impute.bcf.gz"

# Filter alleles to get biallelic, < 10% missing and >= 0.05 MAF and index.
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --include 'F_MISSING < 0.1 && MAF >= 0.05' \
        --min-alleles 2 --max-alleles 2 \
        "$merge_full_bcf" > "$merge_filt_bcf"

bcftools \
    index \
        --force \
        "$merge_filt_bcf"

# Get genotype info from the filtered BCF.
bcftools \
    +smpl-stats \
    "$merge_filt_bcf" > "$MERGE_DIR/gt_info.txt"

# Calculate number of variant sites overall.
nrecords=$(bcftools index --nrecords "$merge_filt_bcf".csi)

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
        "$merge_filt_bcf" | \
bcftools \
    reheader \
        --threads 30 \
        --samples-file "$MERGE_DIR/good_samples.txt" | \
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --min-alleles 2 --max-alleles 2 \
        # Keep only variants with <= 5% missing and >= 5% MAF after subset.
        --include 'F_MISSING <= 0.05 && MAF >= 0.05' \
        > "$merge_filtsub_bcf"

# Get only biallelic SNPs from the biallelic variants bcf.
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --include 'TYPE="snp" && strlen(REF)==1 && strlen(ALT)==1' \
        "$merge_filtsub_bcf" > "$merge_basnpsmaf_bcf"

# Convert to VCF format for imputation with beagle.
bcftools \
    view \
        --threads 30 --no-version --output-type z \
        "$merge_basnpsmaf_bcf" > "$merge_basnpsmaf_vcf"

# Impute missing alleles with Beagle rather than using mean imputation.
# Assumes "beagle" is an alias for java -jar /path/to/beagle.<version>.jar
beagle \
    gt="$merge_basnpsmaf_vcf" \
    out="$merge_basnpsmaf_impute_vcf" \
    nthreads=28

# Re-header the output and convert to bcf.
bcftools \
    view \
        --header-only \
            "$merge_basnpsmaf_bcf" | \
grep "fileformat" > imputed_header.txt

bcftools \
    view \
        --header-only \
            "$merge_basnpsmaf_impute_vcf".vcf.gz | \
grep -P "##INFO|##FORMAT" | grep -v "GT" >> imputed_header.txt

bcftools \
    view \
        --header-only \
            "$merge_basnpsmaf_bcf" | awk 'NR > 1' \
            >> imputed_header.txt

bcftools \
    reheader \
        --header imputed_header.txt \
        "$merge_basnpsmaf_impute_vcf".vcf.gz | \
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        > "$merge_basnpsmaf_impute_bcf"

```

## Estimation of kinship
We use the variants from the previous steps to estimate kinship. We assume the pedigree founder IDs are available in `data/pedigree/inputs/pedig_founder_ids.txt`. Also assumes the file and directory variables set above are still available (e.g. `$merge_basnpsmaf_impute_bcf` etc.).
```bash

# Set $SOFTWARE to local siftware install directory.
export PATH="$SOFTWARE/plink/plink2":"$SOFTWARE/bcftools/plugins/bcftools/":"$SOFTWARE/gcta/gcta-1.95.1-linux-x86_64":"$PATH"
export BCFTOOLS_PLUGINS="$SOFTWARE/bcftools/plugins/bcftools/plugins/"

# Data directory variables.
PEDIGREE_IN_DIR="$WF_DIR/data/pedigree/inputs/"
PEDIGREE_OUT_DIR="$WF_DIR/data/pedigree/out/"
mkdir -p "$PEDIGREE_OUT_DIR"

# Additional input / output file variables.
pedig_founder_ids="$PEDIGREE_IN_DIR/pedig_founder_ids.txt" # (Created with awk '$NF == "founder"' doc/info/varieties_key_pedigree.txt | cut -f 1 > data/pedigree/inputs/pedig_founder_ids.txt)
pedig_cultivar_ids="$PEDIGREE_IN_DIR/pedig_cultivar_ids.txt" # (Created with awk '$NF == "cultivar"' doc/info/varieties_key_pedigree.txt | cut -f 1 > data/pedigree/inputs/pedig_cultivar_ids.txt)
founder_basnpsmaf_impute_bcf="$PEDIGREE_IN_DIR/founder_basnpsmaf_impute.bcf.gz"
pedig_all_ids="$PEDIGREE_IN_DIR/pedig_all_ids.txt"
founder_pos="$PEDIGREE_IN_DIR/founder_pos.txt"
pedig_ldfilt_bcf="$PEDIGREE_IN_DIR/pedig_ldfilt.bcf.gz"
pedig_pos="$PEDIGREE_IN_DIR/pedig_pos.txt"
founders_af="$PEDIGREE_IN_DIR/founders_af.txt"
founders_winter_af="$PEDIGREE_IN_DIR/founders_winter_af.txt"
founders_semiwinter_af="$PEDIGREE_IN_DIR/founders_semiwinter_af.txt"
founders_spring_af="$PEDIGREE_IN_DIR/founders_spring_af.txt"
founders_all_af=$PEDIGREE_IN_DIR/founders_all_af.txt
pedig_base_af="$PEDIGREE_IN_DIR/pedig_base_af.txt"
pedig_ldfilt_plink="$PEDIGREE_IN_DIR/pedig_ldfilt_plink"

# Extract variants biallelic in founders.
bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --samples-file "$pedig_founder_ids" \
        --trim-alt-alleles \
        --output-type u \
        "$merge_basnpsmaf_impute_bcf" | \
bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --max-alleles 2 \
        --min-alleles 2 \
        > "$founder_basnpsmaf_impute_bcf"

cat "$pedig_founder_ids" "$pedig_cultivar_ids" > "$pedig_all_ids"

# Extract these sites from pedigree lines and filter based on linkage disequilibrium.
bcftools \
    query \
        --format '%CHROM\t%POS' \
        "$founder_basnpsmaf_impute_bcf" > "$founder_pos"

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --samples-file "$pedig_all_ids" \
        --targets-file "$founder_pos" \
        --trim-alt-alleles \
        --output-type u \
        "$merge_basnpsmaf_impute_bcf" | \
bcftools \
    +prune \
        --no-version --output-type b \
        --max r2=0.6 \
        > "$pedig_ldfilt_bcf"


# Get these sites from the founders and calculate allele frequencies.
bcftools \
    query \
        --format '%CHROM\t%POS' \
        "$pedig_ldfilt_bcf" > "$pedig_pos"

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --targets-file "$pedig_pos" \
        "$founder_basnpsmaf_impute_bcf" | \
bcftools \
    +fill-tags \
        --output-type u \
        -- --tags AF | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%AF\n' \
        > "$founders_af"

# Subset founders based on their types in the info file.
awk '$1 != "NA" && $3 == "winter"' "$WF_DIR/doc/info/founder_types.txt" | cut -f 1 > "$SUBSET_DIR/founders_winter.txt"
awk '$1 != "NA" && $3 == "semiwinter"' "$WF_DIR/doc/info/founder_types.txt" | cut -f 1 > "$SUBSET_DIR/founders_semiwinter.txt"
awk '$1 != "NA" && $3 == "spring"' "$WF_DIR/doc/info/founder_types.txt" | cut -f 1 > "$SUBSET_DIR/founders_spring.txt"

# Get AFs for base populations.
bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --targets-file "$pedig_pos" \
        --samples-file "$SUBSET_DIR/founders_winter.txt" \
        "$founder_basnpsmaf_impute_bcf" | \
bcftools \
    +fill-tags \
        --output-type u \
        -- --tags AF | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%AF\n' \
        > "$founders_winter_af"

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --targets-file "$pedig_pos" \
        --samples-file "$SUBSET_DIR/founders_semiwinter.txt" \
        "$founder_basnpsmaf_impute_bcf" | \
bcftools \
    +fill-tags \
        --output-type u \
        -- --tags AF | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%AF\n' \
        > "$founders_semiwinter_af"

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --targets-file "$pedig_pos" \
        --samples-file "$SUBSET_DIR/founders_spring.txt" \
        "$founder_basnpsmaf_impute_bcf" | \
bcftools \
    +fill-tags \
        --output-type u \
        -- --tags AF | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%AF\n' \
        > "$founders_spring_af"

# Multiply founder allele frequencies by total base population membership.
winter_af=$(awk 'BEGIN {printf "%.10f", 3/22}')
spring_af=$(awk 'BEGIN {printf "%.10f", 12/22}')
semiwinter_af=$(awk 'BEGIN {printf "%.10f", 7/22}')

# Join the allele frequency files.
paste "$founders_winter_af" "$founders_spring_af" "$founders_semiwinter_af" | \
    awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$9}' > "$founders_all_af"

awk \
    -v waf="$winter_af" -v saf="$spring_af" -v swaf="$semiwinter_af" \
    '{af=$3*waf+$4*saf+$5*swaf}
    af <= 1 {print $1"\t"$2"\t"af}
    af > 1 {print $1"\t"$2"\t"1}' "$founders_all_af" > "$pedig_base_af"

# This makes SNP matrices of just the pedigree lines.
plink2 \
    --bcf "$pedig_ldfilt_bcf" \
    --allow-extra-chr \
    --export A \
    --out "$pedig_ldfilt_plink"

# Get the IDs of the variants for matching allele freqs to SNPs in R.
bcftools \
    view \
        --targets-file "$pedig_pos" \
        "$pedig_ldfilt_bcf" | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%ID' \
        > "$PEDIGREE_IN_DIR/variant_ids.txt"

```
# Calculate VanRaden style kinship
This is the R code used to calculate kinship with the base allele frequencies assumed from above.
```R
library(dplyr)
library(data.table)

# require(tidyr)
# require(ggplot2)
# require(ggrepel)
# require(pheatmap)
# require(cowplot)
# require(GenomicRanges)
# require(vcfR)
# require(reshape2)

# Read the ld filtered plink files into a SNP matrix.
pedig_ldfilt_plink <- "data/pedigree/inputs/pedig_ldfilt_plink"
snpmat_plink <-
  fread(paste0(pedig_ldfilt_plink, ".raw"), data.table = FALSE)

# Remove non-genotype columns and set row names to sample IDs.
snpmat_extract <- as.matrix(snpmat_plink[, 7:ncol(snpmat_plink)])
rownames(snpmat_extract) <- snpmat_plink$IID

# Count alternate alleles instead of reference alleles.
snpmat_alt <- 2 - snpmat_extract

# Read variant IDs file for variants extracted above.
variant_ids_file <- "data/pedigree/inputs/variant_ids.txt"
variant_ids <-
  read.table(
    variant_ids_file, col.names = c("CHROM", "POS", "ID")
  )

# Read the pedigree base allele frequencies.
pedig_base_af_file <- "data/pedigree/inputs/pedig_base_af.txt"
pedig_base_af <-
  read.table(
    pedig_base_af_file, col.names = c("CHROM", "POS", "AF")
  ) |>
  left_join(variant_ids, by = c("CHROM", "POS"))

# plink2 names columns as CHROM:POS_REF_ALT_counted-allele.
snp_keys <- sub("_.*", "", colnames(snpmat_alt))
colnames(snpmat_alt) <- snp_keys

# Align AF vector to SNP matrix column order.
af_matches <-
  pedig_base_af$AF[
    match(snp_keys, pedig_base_af$ID)
  ]

# Drop any SNPs with no AF match.
keep <- !is.na(af_matches)
snpmat_alt <- snpmat_alt[, keep]
pvect <- af_matches[keep]

# Compute VanRaden GRM manually.
zmat <- sweep(snpmat_alt, 2, 2 * pvect, "-")
denom <- 2 * sum(pvect * (1 - pvect))
kinmat_founder <- tcrossprod(zmat) / denom

# The key to get variety names.
varieties_key_pedigree <-
  read.table(
    "doc/info/varieties_key_pedigree.txt",
    sep = "\t",
    header = TRUE
  )

# Variety IDs - maps the bcf / kinship IDs to the line names in the AGG.
varieties_key <-
  read.csv(
    "doc/info/varieties_key.txt",
    header = TRUE, sep = "\t"
  )

varieties_key <-
  varieties_key |> bind_rows(
    data.frame(
      ID = c("Liho", "ganganF73_ref", "no2127_ref", "quintaA_ref", "shengli3_ref", "zheyou73_ref", "zs11_ref"),
      LINE_NAME = c("Liho", "GanganF73", "No2127", "QuintaA", "Shengli3", "Zheyou73", "ZS11")
    )
  )

# Calculate the sum of the off-diagonal elements of the GRM.
pedigcultiv_ids <-
  varieties_key_pedigree |>
  filter(TYPE == "cultivar") |>
  pull(ID)

kinmat_cultiv <-
  kinmat_founder[pedigcultiv_ids, pedigcultiv_ids]

total_sum <- sum(kinmat_cultiv)
diag_sum <- sum(diag(kinmat_cultiv))
sumg_offdiag <- total_sum - diag_sum
print(sumg_offdiag)
print(nrow(kinmat_cultiv))

# Prints 317.8384 for sumg_offdiag and 34 for number of cultivars.

```
# Fixing the allele frequency selection scan
The problem with the previous test was that the number of valid pairs changes depending on missing data (although all genotypes have a minimum of 95% call rate). This is a revised version of the function to perform the selection scan with more detailed comments.
```bash

cultivpedig_ba_bcf="$PEDIGREE_IN_DIR/cultivpedig_ba.bcf.gz"
founders_ba_bcf="$PEDIGREE_IN_DIR/founders_ba.bcf.gz"
founders_winter_af_full="$PEDIGREE_IN_DIR/founders_winter_af_full.txt"
founders_semiwinter_af_full="$PEDIGREE_IN_DIR/founders_semiwinter_af_full.txt"
founders_spring_af_full="$PEDIGREE_IN_DIR/founders_spring_af_full.txt"
founders_all_af_full="$PEDIGREE_IN_DIR/founders_all_af_full.txt"
pedig_base_af_full="$PEDIGREE_IN_DIR/pedig_base_af_full.txt"
selscan_deltap="$PEDIGREE_IN_DIR/selscan_deltap"

# We need variants for only cultivars in the pedigree.
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --samples-file "$pedig_cultivar_ids" \
        --force-samples \
        "$merge_basnpsmaf_bcf" > "$cultivpedig_ba_bcf"

# Then, get the same variants from just the founder lines.
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --samples-file "$pedig_founder_ids" \
        --force-samples \
        "$merge_basnpsmaf_bcf" > "$founders_ba_bcf"

# Get the positions.

# Get AFs for base populations.
bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --samples-file "$SUBSET_DIR/founders_winter.txt" \
        "$founders_ba_bcf" | \
bcftools \
    +fill-tags \
        --output-type u \
        -- --tags AF | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%AF\n' \
        > "$founders_winter_af_full"

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --samples-file "$SUBSET_DIR/founders_semiwinter.txt" \
        "$founders_ba_bcf" | \
bcftools \
    +fill-tags \
        --output-type u \
        -- --tags AF | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%AF\n' \
        > "$founders_semiwinter_af_full"

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --samples-file "$SUBSET_DIR/founders_spring.txt" \
        "$founders_ba_bcf" | \
bcftools \
    +fill-tags \
        --output-type u \
        -- --tags AF | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%AF\n' \
        > "$founders_spring_af_full"

# Multiply founder allele frequencies by total base population membership.
winter_af=$(awk 'BEGIN {printf "%.10f", 3/22}')
spring_af=$(awk 'BEGIN {printf "%.10f", 12/22}')
semiwinter_af=$(awk 'BEGIN {printf "%.10f", 7/22}')

paste "$founders_winter_af_full" \
    "$founders_spring_af_full" "$founders_semiwinter_af_full" | \
    awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$9}' > "$founders_all_af_full"

awk \
    -v waf="$winter_af" -v saf="$spring_af" -v swaf="$semiwinter_af" \
    '{af=$3*waf+$4*saf+$5*swaf}
    af <= 1 {print $1"\t"$2"\t"af}
    af > 1 {print $1"\t"$2"\t"1}' "$founders_all_af_full" > "$pedig_base_af_full"

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
        n_total = 34 # Number of individuals in the kinship matrix.
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

# Then stream the variants in vcf format to the function.
# Run in parallel as it is quite slow.
chromosomes=($(bcftools query --format '%CHROM\n' "$cultivpedig_ba_bcf" | uniq))

printf "%s\n" "${chromosomes[@]}" | \
    parallel --jobs 19 "grep -w {} $pedig_base_af_full > {}_base_af.txt"

printf "%s\n" "${chromosomes[@]}" | \
    parallel --jobs 19 \
    "bcftools query --format '%CHROM\t%POS[\t%GT]\n' $cultivpedig_ba_bcf > {}_varpos.txt"

printf "%s\n" "${chromosomes[@]}" | \
    parallel --jobs 19 \
    "normalise_deltap {}_base_af.txt {}_varpos.txt > {}_selscan_af.txt"

# Join the output files.
cat *selscan_af.txt | \
awk 'BEGIN{print "CHROM\tPOS\tLOCAL_SHARING\tEXPECTED_SHARING\tDEVIATION\tVALID_N\tP_BASE\tP_MODERN"} $1 != "CHROM" {print $0}' > $PEDIGREE_OUT_DIR/cultivpedig_ba_selscan_af.txt

```

## Interpretation of iHS test
This is the code that was used to run the iHS test.

```bash

export PATH=$PATH:"$SOFTWARE/selscan/selscan/bin/linux/"

pocultivar_ids=$PEDIGREE_IN_DIR/pocultivar_ids.txt
pocultivar_basnpsmaf_impute_vcf=$PEDIGREE_IN_DIR/pocultivar_ldfiltsel.vcf.gz
pocultivar_basnpsmaf_impute_bcf=$PEDIGREE_IN_DIR/pocultivar_basnpsmaf_impute.bcf.gz
pedig_all_basnpsmaf_impute_bcf=$PEDIGREE_IN_DIR/pedig_all_basnpsmaf_impute.bcf.gz

# This is added - run the iHS scan and other diversity tests also on
# Australian breeding lines.
ausbl_ids=$PEDIGREE_IN_DIR/ausbl_ids.txt
ausbl_basnpsmaf_impute_vcf=$PEDIGREE_IN_DIR/ausbl_basnpsmaf_impute.vcf.gz
ausbl_basnpsmaf_impute_bcf=$PEDIGREE_IN_DIR/ausbl_basnpsmaf_impute.bcf.gz

# Get input files with just Australian cultivars and just the rest of the
# cultivars.
awk '$NF ~ /cultivar/ {print $1}' doc/info/varieties_key_pedigree.txt > $pocultivar_ids

bcftools \
    view \
        --no-version --output-type u --threads 30 \
        --samples-file $pocultivar_ids \
        --trim-alt-alleles \
        --force-samples \
        $merge_basnpsmaf_impute_bcf | \
bcftools \
    view \
        --no-version --output-type z --threads 30 \
        --min-alleles 2 --max-alleles 2 \
        > $pocultivar_basnpsmaf_impute_vcf

bcftools \
    view \
        --no-version --output-type u --threads 30 \
        --samples-file $pocultivar_ids \
        --trim-alt-alleles \
        --force-samples \
        $merge_basnpsmaf_impute_bcf | \
bcftools \
    view \
        --no-version --output-type b --threads 30 \
        --min-alleles 2 --max-alleles 2 \
        > $pocultivar_basnpsmaf_impute_bcf

bcftools \
    view \
        --no-version --output-type u --threads 30 \
        --samples-file $pedig_all_ids \
        --trim-alt-alleles \
        --force-samples \
        $merge_basnpsmaf_impute_bcf | \
bcftools \
    view \
        --no-version --output-type b --threads 30 \
        --min-alleles 2 --max-alleles 2 \
        > $pedig_all_basnpsmaf_impute_bcf

bcftools \
    index \
        --force \
        $pocultivar_basnpsmaf_impute_vcf

bcftools \
    index \
        --force \
        $pocultivar_basnpsmaf_impute_bcf

bcftools \
    index \
        --force \
        $pedig_all_basnpsmaf_impute_bcf

# Split input into chromosomes and make dummy maps.
chrs=(A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 C01 C02 C03 C04 C05 C06 C07 C08 C09)

for chr in "${chrs[@]}"; do

    chr_vcf=$PEDIGREE_IN_DIR/pocultivar_basnpsmaf_impute_${chr}.vcf.gz
    chr_bcf=$PEDIGREE_IN_DIR/pocultivar_basnpsmaf_impute_${chr}.bcf.gz
    chr_map=$PEDIGREE_IN_DIR/pocultivar_basnpsmaf_impute_${chr}.map

    # Split by chromosome.
    bcftools view \
        --no-version --output-type b --threads 30 \
        --regions $chr \
        $pocultivar_basnpsmaf_impute_vcf > $chr_bcf

    bcftools index $chr_bcf

    # Make a dummy map: id, position (bp), position (pseudo-cM = bp/1e6).
    bcftools query \
        -f '%ID\t%POS\n' $chr_bcf | \
    awk -v chr=$chr '{print chr"\t"$1"\t"$2/1e6"\t"$2}' \
        > $chr_map

    bcftools view \
        --no-version --output-type z --threads 30 \
        --regions $chr \
        $pocultivar_basnpsmaf_impute_vcf > $chr_vcf

done

# Run selscan on the chromosomes.
for chr in "${chrs[@]}"; do

    chr_vcf=$PEDIGREE_IN_DIR/pocultivar_basnpsmaf_impute_${chr}.vcf.gz
    chr_map=$PEDIGREE_IN_DIR/pocultivar_basnpsmaf_impute_${chr}.map

    selscan \
        --ihs \
        --vcf $chr_vcf \
        --map $chr_map \
        --keep-low-freq \
        --threads 28 \
        --out $PEDIGREE_IN_DIR/ihs_cultiv_${chr}

done

# Normalise across all chromosomes.
selscan \
    norm \
        --ihs \
        --files $(ls $PEDIGREE_IN_DIR/ihs_cultiv_*.ihs.out) \
        --bp-win \
        --winsize 100000

```

# Other tests
```bash
cultivall_ba_bcf=$PEDIGREE_IN_DIR/cultivall_ba.bcf.gz
pedig_ba_bcf=$PEDIGREE_IN_DIR/pedig_ba.bcf.gz
cultivall_ba_pi="$PEDIGREE_OUT_DIR/cultivall_ba_pi.txt"
cultivall_ba_tajima="$PEDIGREE_IN_DIR/cultivall_ba_tajima"

# Subset the filtered bcf.
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --samples-file $pocultivar_ids \
        --force-samples \
        $merge_basnpsmaf_bcf > $cultivall_ba_bcf

bcftools \
    index \
        --force \
        $cultivall_ba_bcf

bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --samples-file $pedig_all_ids \
        --force-samples \
        $merge_basnpsmaf_bcf > $pedig_ba_bcf

bcftools \
    index \
        --force \
        $pedig_ba_bcf

# Nucleotide diversity in sliding windows
vcftools \
    --bcf $cultivall_ba_bcf \
    --window-pi 100000 \
    --window-pi-step 50000 \
    --out $cultivall_ba_pi

# Calculate Tajima's D.
vcftools \
    --bcf $cultivall_ba_bcf \
    --TajimaD 100000 \
    --out $cultivall_ba_tajima

# Input for Nei's haplotype diversity.
acultivar_ids=$PEDIGREE_IN_DIR/acultivar_ids.txt
acultivar_basnpsmaf_bcf=$PEDIGREE_IN_DIR/acultivar_basnpsmaf.bcf.gz

awk 'NR > 1 {print $1}' doc/info/varieties_key_pedigree.txt > $acultivar_ids

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --force-samples --trim-alt-alleles \
        --samples-file $acultivar_ids \
        $merge_basnpsmaf_bcf | \
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --min-alleles 2 \
        > $acultivar_basnpsmaf_bcf

bcftools index --force $acultivar_basnpsmaf_bcf

```
# Assessing the A08 LEA allele transmission through the pedigree
```bash
source envs/bash/revisions.sh

awk '$NF != "cultivar-missing" && NR > 1' doc/info/varieties_key_pedigree.txt \
    > doc/info/pedigree_samples.txt

bcftools index --force $merge_basnpsmaf_bcf

bcftools \
    view \
    --samples-file doc/info/pedigree_samples.txt \
    --force-samples \
    --regions A08 \
    --output-type v \
    $merge_basnpsmaf_bcf

```
