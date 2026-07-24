# Environment setup
```bash
source envs/bash/gwas.sh
```
# Generate input data for BLINK
## Prepare variant matrix
```bash

# Start by filtering variants to remove all perfectly linked alleles.
bcftools \
    +prune \
        --max r2=1.0 \
        --window 1000kb \
        --output-type b \
        $merge_basnpsmaf_impute_bcf \
        > $merge_basnpsmaf_iprune_bcf

# Get matrix of allele counts.
plink2 \
    --bcf $merge_basnpsmaf_iprune_bcf \
    --allow-extra-chr \
    --export A \
    --make-bed \
    --out $merge_basnpsmaf_iprune_plink

# Reformat for blink.
cut -f 7- $merge_basnpsmaf_iprune_plink.raw | \
awk 'NR > 1' > $merge_basnpsmaf_iprune_snpmat

sed -e 's/^A01\t/1\t/' -e 's/^A02\t/2\t/' -e 's/^A03\t/3\t/' \
    -e 's/^A04\t/4\t/' -e 's/^A05\t/5\t/' -e 's/^A06\t/6\t/' \
    -e 's/^A07\t/7\t/' -e 's/^A08\t/8\t/' -e 's/^A09\t/9\t/' \
    -e 's/^A10\t/10\t/' -e 's/^C01\t/11\t/' -e 's/^C02\t/12\t/' \
    -e 's/^C03\t/13\t/' -e 's/^C04\t/14\t/' -e 's/^C05\t/15\t/' \
    -e 's/^C06\t/16\t/' -e 's/^C07\t/17\t/' -e 's/^C08\t/18\t/' \
    -e 's/^C09\t/19\t/' $merge_basnpsmaf_iprune_plink.bim | \
awk 'BEGIN {print "rs\tchr\tpos"} {print $2"\t"$1"\t"$4}' \
    > $merge_basnpsmaf_iprune_map

cut -f 2 $merge_basnpsmaf_iprune_plink.raw | awk 'NR > 1' \
    > $merge_basnpsmaf_iprune_ids

```
## Prepare relationship matrix for PCA
```bash

# Filter a second time for LD to get relationship matrices for
# PCA.
bcftools \
    +prune \
        --max r2=1.0 \
        --window 1000kb \
        --output-type b \
        $merge_basnpsmaf_iprune_bcf \
        > $merge_basnpsmaf_ldfilt_bcf

# Get matrix of allele counts.
plink2 \
    --bcf $merge_basnpsmaf_ldfilt_bcf \
    --allow-extra-chr \
    --export A \
    --make-bed \
    --out $merge_basnpsmaf_ldfilt_plink

awk 'BEGIN {
    split("A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 C01 C02 C03 C04 C05 C06 C07 C08 C09", chrs, " ")
    for (i=1; i<=length(chrs); i++) map[chrs[i]] = i
}
{ $1 = map[$1]; print }' OFS='\t' $merge_basnpsmaf_ldfilt_plink.bim \
    > tmp.bim && mv tmp.bim $merge_basnpsmaf_ldfilt_plink.bim

gcta64 \
    --bfile $merge_basnpsmaf_ldfilt_plink \
    --make-grm-alg 1 \
    --make-grm \
    --out $merge_basnpsmaf_ldfilt_grm

```