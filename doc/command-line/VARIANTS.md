# Dependencies
vg 1.70.0
bgzip 1.11
bcftools 1.9
htslib 1.9

# Code

## Set up environment

```bash
source envs/bash/variants.sh
```

## Call and filter variants directly from the graph

Using vg deconstruct, get the variants from the graph.

```bash

# Call directly from the graph.
vg \
    deconstruct \
        --all-snarls \
        --path-prefix $ref_name \
        --snarls $graph_snarls \
        --threads 30 \
        $graph_gbz > $graph_vcf

# Set conflicting genotypes to missing, rename
# chromosomes, remove INFO except AT and all
# FILTER IDs.
./util/reannotate_vcf_conflicts.py \
    $graph_vcf | \
bcftools \
    annotate \
        --threads 30 --no-version --output-type b \
        --include 'LV=0' --remove \
            INFO/AC,INFO/AF,INFO/NS,INFO/AN,FILTER \
        --rename-chrs $chromids_pansn_file \
        - > $graph_bcf

```

This variant file contains 8,014,296 sites.

## Call and filter variants from reads mapped to the graph

We take the merged vcf from the graphcall workflow 'workflows/graphcall.nf' and filter it.
Heterozygous genotypes and genotypes with GQ < 20 or DP < 2 are set to missing.

```bash

# Set heterozygous, low quality (GQ < 20 + DP < 2) genotypes to missing
# and rename chromosomes for reads.
bcftools \
    +setGT \
        --threads 30 --no-version --output-type u \
        $reads_vcf \
        -- \
        --include 'GT == "het" | FMT/GQ < 20 | FMT/DP < 2' \
        --target-gt q --new-gt './.' | \
bcftools \
    annotate \
        --threads 30 --no-version --output-type b \
        --remove INFO/LV,INFO/PS,INFO/DP,FILTER,FORMAT/AD,FORMAT/MAD,FORMAT/DP,FORMAT/GL,FORMAT/GQ,FORMAT/GP,FORMAT/XD --rename-chrs $chromids_file \
        - > $reads_filt_bcf

```

## Merge variants from direct graph calling and read mapping based calling

The two files are simply merged, after removing all duplicates. Alleles differing between reads and graph are ignored. The final filtering procedure removes variants with missing data in >= 10% of lines and a MAF < 0.01.

```bash

# Merge graph and reads bcfs.
bcftools \
    index \
        --threads 30 $graph_bcf

bcftools \
    index \
        --threads 30 $reads_filt_bcf

# Initial naive merge.
bcftools \
    merge \
        --threads 30 --no-version --output-type u \
        --merge none \
        $graph_bcf $reads_filt_bcf | \
bcftools \
    norm \
        --threads 30 --no-version --output-type b \
        --rm-dup all \
        - > $merge_full_bcf

```

## Fix chromosome IDs
I made a mistake converting X5A accessions to chromosome names. I have used the minigraph cactus alignments<\br>
to determine which chromosome is which. This is run once, then all subsequent files are re-created.
```bash

mkdir -p tmp

# Get a file for renaming chromosomes.
awk 'NR > 1 {print $3"\t"$2}'  doc/info/chromosomes_rename_key.txt > tmp/chrom_rename.txt

# Create a temporary file with renamed chromosomes
# then convert it back to the full bcf.
bcftools \
    annotate \
        --threads 30 --no-version --output-type b \
        --rename-chrs tmp/chrom_rename.txt \
        $merge_full_bcf > tmp.bcf.gz

mv tmp.bcf.gz $merge_full_bcf

```

## Filter variants and samples based on missingness and MAF

``` bash

# Remove variants with too much missing data and a low MAF.
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --include 'F_MISSING < 0.1 && MAF >= 0.05' \
        --min-alleles 2 --max-alleles 2 \
        $merge_full_bcf > $merge_filt_bcf

# Get information on each sample.
bcftools \
    +smpl-stats \
        $merge_filt_bcf > $MERGE_DIR/gt_info.txt

bcftools \
    index \
        --force \
        $merge_filt_bcf

# Calculate number of variant sites overall.
nrecords=$(bcftools index --nrecords $merge_filt_bcf.csi)

# Get missing rate.
awk '$1 == "FLT0"' $MERGE_DIR/gt_info.txt | \
    cut -f 2,12 | awk -v nrec=$nrecords '{print $1"\t"$2/nrec}' \
    > $MERGE_DIR/gt_missrate.txt

awk '$2 <= 0.26' $MERGE_DIR/gt_missrate.txt | \
    cut -f 1 > $MERGE_DIR/good_samples.txt

# Subset to better quality samples.
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
        --include 'F_MISSING <= 0.05 && MAF >= 0.05' \
        > $merge_filtsub_bcf

# Subset to better quality samples - no MAF filtering here. Mainly used
# for targeted rare allele analyses.
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
        --include 'F_MISSING <= 0.05' \
        > $merge_ba_bcf

# Get biallelic SNPs from the biallelic variants bcf.
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --include 'TYPE="snp" && strlen(REF)==1 && strlen(ALT)==1' \
        $merge_filtsub_bcf > $merge_basnpsmaf_bcf

# Then SVs.
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --include 'strlen(REF)>=50 || strlen(ALT)>=50' \
        $merge_filtsub_bcf > $merge_basvs_bcf
```

## Imputation of missing alleles

```bash

# Convert biallelic filtered SNPs to vcf.
bcftools \
    view \
        --threads 30 --no-version --output-type z \
        $merge_basnpsmaf_bcf > $merge_basnpsmaf_vcf

# Impute missing alleles with Beagle rather than using mean imputation.
beagle \
    gt=$merge_basnpsmaf_vcf \
    out=$merge_basnpsmaf_impute_vcf \
    nthreads=28

# Re-header the output and convert to bcf.
bcftools \
    view \
        --header-only \
            $merge_basnpsmaf_bcf | \
grep "fileformat" > imputed_header.txt

bcftools \
    view \
        --header-only \
            $merge_basnpsmaf_impute_vcf.vcf.gz | \
grep -P "##INFO|##FORMAT" | grep -v "GT" >> imputed_header.txt

bcftools \
    view \
        --header-only \
            $merge_basnpsmaf_bcf | awk 'NR > 1' \
            >> imputed_header.txt

bcftools \
    reheader \
        --header imputed_header.txt \
        $merge_basnpsmaf_impute_vcf.vcf.gz | \
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        > $merge_basnpsmaf_impute_bcf

```

## LD pruning for genomic prediction and PCA

```bash

# Filter to remove snps that are in 0.6 r2 LD for PCA.
bcftools \
    +prune \
        --no-version --output-type b \
        --max r2=0.6 \
        $merge_basnpsmaf_impute_bcf > $merge_ldfiltpca_bcf

# Filter to remove SNPs in very high LD 0.99 for computational efficiency of GP.
bcftools \
    +prune \
        --no-version --output-type b \
        --window 50kb \
        --max r2=0.99 \
        $merge_basnpsmaf_impute_bcf > $merge_ldfiltgp_bcf

```

## Run PCA for population structure covariates

```bash

# Convert the LD-filtered BCF to plink format.
plink2 \
    --bcf $merge_ldfiltpca_bcf \
    --allow-extra-chr \
    --make-bed \
    --out $merge_ldfiltpca_plink

awk 'BEGIN {
    split("A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 C01 C02 C03 C04 C05 C06 C07 C08 C09", chrs, " ")
    for (i=1; i<=length(chrs); i++) map[chrs[i]] = i
}
{ $1 = map[$1]; print }' OFS='\t' $merge_ldfiltpca_plink.bim > tmp.bim && mv tmp.bim $merge_ldfiltpca_plink.bim

# Create the vanraden GRM with gcta.
gcta64 \
    --bfile $merge_ldfiltpca_plink \
    --make-grm-alg 1 \
    --make-grm \
    --out $merge_ldfiltpca_grm

```
## Make input matrices for genomic prediction

```bash

# Convert the LD-filtered GP BCF to plink format.
plink2 \
    --bcf $merge_ldfiltgp_bcf \
    --allow-extra-chr \
    --make-bed \
    --out $merge_ldfiltgp_plink

awk 'BEGIN {
    split("A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 C01 C02 C03 C04 C05 C06 C07 C08 C09", chrs, " ")
    for (i=1; i<=length(chrs); i++) map[chrs[i]] = i
}
{ $1 = map[$1]; print }' OFS='\t' $merge_ldfiltgp_plink.bim > tmp.bim && mv tmp.bim $merge_ldfiltgp_plink.bim

# Create the vanraden GRM with gcta.
gcta64 \
    --bfile $merge_ldfiltgp_plink \
    --make-grm \
    --make-grm-alg 1 \
    --out $merge_ldfiltgp_grm

```

# Prepare input for BLINK.
# Impute missing alleles with beagle.
bcftools \
    filter \
    --include 'N_ALT == 1' \
    --output-type u \
    --threads 30 \
    $gwas_bcf_full | \
bcftools \
    norm \
    --threads 30 \
    --rm-dup all \
    --output-type z \
    - > $gwas_vcf_full

beagle \
    gt=$gwas_vcf_full \
    out=$gwas_vcf_impute \
    nthreads=15

bcftools \
    reheader \
    --fai $ref_idx_chrids \
    --threads 30 \
    $gwas_vcf_impute | \
bcftools \
    view \
    --output-type b \
    - > $gwas_bcf_impute.bcf.gz

plink2 \
    --export A \
    --bcf $gwas_bcf_impute.bcf.gz \
    --out $gwas_inputs/gwas_blink_gt

# Subset the imputed bcf to get 60K markers.
# Build GRMs for genomic prediction from filtered
# and non-filtered.
mkdir -p data/gp/inputs/
gblup_kinship_full="data/gp/inputs/gpin_kinship_full"
gblup_kinship_60k="data/gp/inputs/gpin_kinship_60k"

plink2 \
    --bcf $gwas_bcf_full \
    --make-rel square \
    --out $gblup_kinship_full

bcftools \
    query \
        --format "%CHROM\t%POS" \
        $gwas_bcf_full | \
shuf \
    -n 60000 \
    - > data/gp/inputs/random_sites60k.txt

bcftools \
    index \
        --threads 30 \
        --force \
        $gwas_bcf_full

bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --regions-file data/gp/inputs/random_sites60k.txt \
        $gwas_bcf_full > data/gwas/inputs/gwasin_60k.bcf.gz

plink2 \
    --bcf data/gwas/inputs/gwasin_60k.bcf.gz \
    --make-rel square \
    --out $gblup_kinship_60k

# This is for input to BLINK.
cut -f 7- $gwas_inputs/gwas_blink_gt.raw | awk 'NR > 1' \
    > $gwas_inputs/gwas_blink.gt
cut -f 2 $gwas_inputs/gwas_blink_gt.raw | awk 'NR > 1' \
    > $gwas_inputs/sample_ids.txt

bcftools \
    view \
        --output-type v \
        $gwas_bcf_impute.full.bcf.gz | \
grep \
    -v "^#" | \
awk 'BEGIN {print "rs\tchr\tpos"}
           {print $3"\t"$1"\t"$2}' \
    > $gwas_inputs/gwas_blink.map

```

