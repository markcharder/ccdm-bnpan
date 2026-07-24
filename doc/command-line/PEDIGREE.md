# Dependencies
bcftools 1.9
htslib 1.9

# Environment
```bash
source envs/bash/pedigree.sh
```

# Code
## Estimating diversity loss from founders to cultivars
Here, I aim to assess formally how much diversity is lost from founders to cultivars. My assumption is that the founders<\br>
we have represent less diversity than there was in the full set, since the missing lines are known to be divergent from<\br>
the ones we do have. For example, two are B. juncea.<\br>
For the analysis, I will split the groups up into founders, pedigree cultivars and non-pedigree cultivars. Keeping the latter<\br>
as a second group will help assess whether diversity is higher or lower in new pedigree-adjacent cultivars than pedigree<\br>
cultivars. If it is substantially higher on average, it is possible that new alleles have been brought in.<\br>
The first step is to extract founders, pedigree cultivars and non-pedigree cultivars from the filtered biallelic SNPs file.<\br>
```bash

awk \
    -v pfids=$pedig_founder_ids \
    -v pcids=$pedig_cultivar_ids \
    -v ocids=$other_cultivar_ids \
    'NR > 1 && $NF == "founder" {print $1 > pfids}
     NR > 1 && $NF == "cultivar" {print $1 > pcids}
     NR > 1 && $NF == "cultivar-missing" {print $1 > ocids}' doc/info/varieties_key_pedigree.txt

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --force-samples --trim-alt-alleles \
        --samples-file $pedig_founder_ids \
        $merge_basnpsmaf_bcf | \
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --min-alleles 2 \
        > $pfounder_basnpsmaf_bcf

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --force-samples --trim-alt-alleles \
        --samples-file $pedig_cultivar_ids \
        $merge_basnpsmaf_bcf | \
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --min-alleles 2 \
        > $pcultivar_basnpsmaf_bcf

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --force-samples --trim-alt-alleles \
        --samples-file $pedig_cultivar_ids \
        $merge_basnpsmaf_bcf | \
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --min-alleles 2 \
        > $ocultivar_basnpsmaf_bcf

bcftools index --force $pfounder_basnpsmaf_bcf
bcftools index --force $pcultivar_basnpsmaf_bcf
bcftools index --force $ocultivar_basnpsmaf_bcf

# Get one with all cultivars and founders in.
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

Then, for each chromosome, calculate Watterson's theta in the three subsets.
```bash

n_pfounder=$(bcftools query --list-samples $pfounder_basnpsmaf_bcf | wc -l)
n_pcultivar=$(bcftools query --list-samples $pcultivar_basnpsmaf_bcf | wc -l)
n_ocultivar=$(bcftools query --list-samples $ocultivar_basnpsmaf_bcf | wc -l)

chromlens=($(bcftools \
    view \
        --header-only $merge_basnpsmaf_bcf | \
        grep "^##contig" | \
        sed 's/.*ID=\([^,>]*\).*length=\([0-9]*\).*/\1\t\2/' | \
        sort -V -k1,1 | \
        cut -f 2))

for ((i=0;i<${#chromosomes[@]};i++)); do

    chrom=${chromosomes[$i]}
    clen=${chromlens[$i]}

    s_pfounder=$(bcftools stats --regions $chrom $pfounder_basnpsmaf_bcf | awk '/^SN.*number of SNPs/{print $NF}')
    s_pcultivar=$(bcftools stats --regions $chrom $pcultivar_basnpsmaf_bcf | awk '/^SN.*number of SNPs/{print $NF}')
    s_ocultivar=$(bcftools stats --regions $chrom $ocultivar_basnpsmaf_bcf | awk '/^SN.*number of SNPs/{print $NF}')

    echo -e "$chrom\t$s_pfounder\t$s_pcultivar\t$s_ocultivar\t$n_pfounder\t$n_pcultivar\t$n_ocultivar\t$clen"
done | \
awk 'BEGIN {print "CHROM\tPFOUNDER\tPCULTIVAR\tOCULTIVAR"}
    NR == 1 {
        anpf = 0; for (i = 1; i <= (2*$5-1); i++) {anpf += 1/i}
        anpc = 0 ; for (i = 1; i <= (2*$6-1); i++) {anpc += 1/i}
        anoc = 0; for (i = 1; i <= (2*$7-1); i++) {anoc += 1/i}
    }
    {
        printf "%s\t%.8f\t%.8f\t%.8f\n", $1, $2/($8*anpf), $3/($8*anpc), $4/($8*anoc)
    }' > $PEDIGREE_OUT_DIR/theta_w_pedigree.txt

```
Then calculate loss of segregating sites from founders.
```bash

for ((i=0;i<${#chromosomes[@]};i++)); do
    chrom=${chromosomes[$i]}

    if (($i == 0)); then
        echo -e "CHROM\tSFOUNDER\tSRETAINED"
    fi

    s_founder=$(bcftools view --regions $chrom $pfounder_basnpsmaf_bcf | \
        bcftools view --min-alleles 2 | bcftools stats | awk '/^SN.*number of SNPs/{print $NF}')

    bcftools view --min-alleles 2 --regions $chrom $pfounder_basnpsmaf_bcf | \
        bcftools query --format '%CHROM\t%POS\n' > founder_sites.txt
    s_retained=$(bcftools view --regions $chrom --targets-file founder_sites.txt $pcultivar_basnpsmaf_bcf | \
        bcftools view --min-alleles 2 | bcftools stats | awk '/^SN.*number of SNPs/{print $NF}')

    echo -e "$chrom\t$s_founder\t$s_retained"

done > $PEDIGREE_OUT_DIR/alleles_retained.txt

```
## Estimate structural variant loss across chromosomes
Then, I am going to estimate loss of structural variants from founders to cultivars.
```bash

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --force-samples --trim-alt-alleles \
        --samples-file $pedig_founder_ids \
        $merge_basvs_bcf | \
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --min-alleles 2 \
        > $pfounder_basvs_bcf

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --force-samples --trim-alt-alleles \
        --samples-file $pedig_cultivar_ids \
        $merge_basvs_bcf | \
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --min-alleles 2 \
        > $pcultivar_basvs_bcf

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --force-samples --trim-alt-alleles \
        --samples-file $other_cultivar_ids \
        $merge_basvs_bcf | \
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --min-alleles 2 \
        > $ocultivar_basvs_bcf

bcftools index --force $pfounder_basvs_bcf
bcftools index --force $pcultivar_basvs_bcf
bcftools index --force $ocultivar_basvs_bcf

```

Then, for each chromosome, calculate Watterson's theta in the three subsets.
```bash

n_pfounder=$(bcftools query --list-samples $pfounder_basvs_bcf | wc -l)
n_pcultivar=$(bcftools query --list-samples $pcultivar_basvs_bcf | wc -l)
n_ocultivar=$(bcftools query --list-samples $ocultivar_basvs_bcf | wc -l)

chromlens=($(bcftools \
    view \
        --header-only $merge_basvs_bcf | \
        grep "^##contig" | \
        sed 's/.*ID=\([^,>]*\).*length=\([0-9]*\).*/\1\t\2/' | \
        sort -V -k1,1 | \
        cut -f 2))

for ((i=0;i<${#chromosomes[@]};i++)); do

    chrom=${chromosomes[$i]}
    clen=${chromlens[$i]}

    s_pfounder=$(bcftools stats --regions $chrom $pfounder_basvs_bcf | awk '/^SN.*number of records/{print $NF}')
    s_pcultivar=$(bcftools stats --regions $chrom $pcultivar_basvs_bcf | awk '/^SN.*number of records/{print $NF}')
    s_ocultivar=$(bcftools stats --regions $chrom $ocultivar_basvs_bcf | awk '/^SN.*number of records/{print $NF}')

    echo -e "$chrom\t$s_pfounder\t$s_pcultivar\t$s_ocultivar\t$n_pfounder\t$n_pcultivar\t$n_ocultivar\t$clen\t$s_pfounder\t$s_pcultivar\t$s_ocultivar"
done | \
awk 'BEGIN {print "CHROM\tPFOUNDER\tPCULTIVAR\tOCULTIVAR"}
    NR == 1 {
        anpf = 0; for (i = 1; i <= (2*$5-1); i++) {anpf += 1/i}
        anpc = 0 ; for (i = 1; i <= (2*$6-1); i++) {anpc += 1/i}
        anoc = 0; for (i = 1; i <= (2*$7-1); i++) {anoc += 1/i}
    }
    {
        printf "%s\t%.10f\t%.10f\t%.10f\n", $1, $2/($8*anpf), $3/($8*anpc), $4/($8*anoc)
    }' > $PEDIGREE_OUT_DIR/theta_w_pedigree_sv.txt

```
Then calculate loss of segregating sites from founders.
```bash

for ((i=0;i<${#chromosomes[@]};i++)); do
    chrom=${chromosomes[$i]}

    if (($i == 0)); then
        echo -e "CHROM\tSFOUNDER\tSRETAINED"
    fi

    s_founder=$(bcftools view --regions $chrom $pfounder_basvs_bcf | \
        bcftools view --min-alleles 2 | bcftools stats | awk '/^SN.*number of records/{print $NF}')

    bcftools view --min-alleles 2 --regions $chrom $pfounder_basvs_bcf | \
        bcftools query --format '%CHROM\t%POS\n' > founder_sites.txt
    s_retained=$(bcftools view --regions $chrom --targets-file founder_sites.txt $pcultivar_basvs_bcf | \
        bcftools view --min-alleles 2 | bcftools stats | awk '/^SN.*number of records/{print $NF}')

    echo -e "$chrom\t$s_founder\t$s_retained"

done > $PEDIGREE_OUT_DIR/alleles_retained_sv.txt
```

## Get sites for estimating relationships between lines in the pedigrees

```bash

# Get all imputed sites biallelic in founders.
bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --samples-file $pedig_founder_ids \
        --trim-alt-alleles \
        --output-type u \
        $merge_basnpsmaf_impute_bcf | \
bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --max-alleles 2 \
        --min-alleles 2 \
        > $founder_basnpsmaf_impute_bcf

cat $pedig_founder_ids $pedig_cultivar_ids > $pedig_all_ids

# Extract these sites from pedigree lines and filter based on linkage disequilibrium.
bcftools \
    query \
        --format '%CHROM\t%POS' \
        $founder_basnpsmaf_impute_bcf > $founder_pos

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --samples-file $pedig_all_ids \
        --targets-file $founder_pos \
        --trim-alt-alleles \
        --output-type u \
        $merge_basnpsmaf_impute_bcf | \
bcftools \
    +prune \
        --no-version --output-type b \
        --max r2=0.6 \
        > $pedig_ldfilt_bcf

# Get these sites from the founders and calculate allele frequencies.
bcftools \
    query \
        --format '%CHROM\t%POS' \
        $pedig_ldfilt_bcf > $pedig_pos

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --targets-file $pedig_pos \
        $founder_basnpsmaf_impute_bcf | \
bcftools \
    +fill-tags \
        --output-type u \
        -- --tags AF | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%AF\n' \
        > $founders_af

# Subset founders based on their types in the info file.
awk '$1 != "NA" && $3 == "winter"' doc/info/founder_types.txt | cut -f 1 > $SUBSET_DIR/founders_winter.txt
awk '$1 != "NA" && $3 == "semiwinter"' doc/info/founder_types.txt | cut -f 1 > $SUBSET_DIR/founders_semiwinter.txt
awk '$1 != "NA" && $3 == "spring"' doc/info/founder_types.txt | cut -f 1 > $SUBSET_DIR/founders_spring.txt

# Get AFs for base populations.
bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --targets-file $pedig_pos \
        --samples-file $SUBSET_DIR/founders_winter.txt \
        $founder_basnpsmaf_impute_bcf | \
bcftools \
    +fill-tags \
        --output-type u \
        -- --tags AF | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%AF\n' \
        > $founders_winter_af

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --targets-file $pedig_pos \
        --samples-file $SUBSET_DIR/founders_semiwinter.txt \
        $founder_basnpsmaf_impute_bcf | \
bcftools \
    +fill-tags \
        --output-type u \
        -- --tags AF | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%AF\n' \
        > $founders_semiwinter_af

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --targets-file $pedig_pos \
        --samples-file $SUBSET_DIR/founders_spring.txt \
        $founder_basnpsmaf_impute_bcf | \
bcftools \
    +fill-tags \
        --output-type u \
        -- --tags AF | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%AF\n' \
        > $founders_spring_af

# Multiply founder allele frequencies by total base population membership.
winter_af=$(awk 'BEGIN {printf "%.10f", 3/22}')
spring_af=$(awk 'BEGIN {printf "%.10f", 12/22}')
semiwinter_af=$(awk 'BEGIN {printf "%.10f", 7/22}')

paste $founders_winter_af $founders_spring_af $founders_semiwinter_af | \
    awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$9}' > $founders_all_af

awk \
    -v waf=$winter_af -v saf=$spring_af -v swaf=$semiwinter_af \
    '{af=$3*waf+$4*saf+$5*swaf}
    af <= 1 {print $1"\t"$2"\t"af}
    af > 1 {print $1"\t"$2"\t"1}' $founders_all_af > $pedig_base_af

# This makes SNP matrices of just the pedigree lines.
plink2 \
    --bcf $pedig_ldfilt_bcf \
    --allow-extra-chr \
    --export A \
    --out $pedig_ldfilt_plink

# Get the IDs of the variants for matching allele freqs to SNPs in R.
bcftools \
    view \
        --targets-file $pedig_pos \
        $pedig_ldfilt_bcf | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%ID' \
        > $PEDIGREE_IN_DIR/variant_ids.txt

```

# Selection scans of Australian cultivars

## iHS test using selscan

```bash

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

## Custom allele frequency change test
```bash

# We need variants for only cultivars in the pedigree.
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --samples-file $pedig_cultivar_ids \
        --force-samples \
        $merge_ba_bcf > $cultivpedig_ba_bcf

# Then, get the same variants from just the founder lines.
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --samples-file $pedig_founder_ids \
        --force-samples \
        $merge_ba_bcf > $founders_ba_bcf

# Get the positions.

# Get AFs for base populations.
bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --samples-file $SUBSET_DIR/founders_winter.txt \
        $founders_ba_bcf | \
bcftools \
    +fill-tags \
        --output-type u \
        -- --tags AF | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%AF\n' \
        > $founders_winter_af_full

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --samples-file $SUBSET_DIR/founders_semiwinter.txt \
        $founders_ba_bcf | \
bcftools \
    +fill-tags \
        --output-type u \
        -- --tags AF | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%AF\n' \
        > $founders_semiwinter_af_full

bcftools \
    view \
        --threads 30 --no-version --output-type u \
        --samples-file $SUBSET_DIR/founders_spring.txt \
        $founders_ba_bcf | \
bcftools \
    +fill-tags \
        --output-type u \
        -- --tags AF | \
bcftools \
    query \
        --format '%CHROM\t%POS\t%AF\n' \
        > $founders_spring_af_full

# Multiply founder allele frequencies by total base population membership.
winter_af=$(awk 'BEGIN {printf "%.10f", 3/22}')
spring_af=$(awk 'BEGIN {printf "%.10f", 12/22}')
semiwinter_af=$(awk 'BEGIN {printf "%.10f", 7/22}')

paste $founders_winter_af_full \
    $founders_spring_af_full $founders_semiwinter_af_full | \
    awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$9}' > $founders_all_af_full

awk \
    -v waf=$winter_af -v saf=$spring_af -v swaf=$semiwinter_af \
    '{af=$3*waf+$4*saf+$5*swaf}
    af <= 1 {print $1"\t"$2"\t"af}
    af > 1 {print $1"\t"$2"\t"1}' $founders_all_af_full > $pedig_base_af_full

# Calculate base allele frequencies as before.

# Function for estimating expected deviation from global relatedness for each SNP.
normalise_deltap() {
    
    base_freqs=$1
    varpos=$2

    awk -v base_freqs=$base_freqs '
    BEGIN {

        FS = "\t";
        OFS = "\t";

        sum_G_offdiag = 317.8384; 
    
        while ((getline < base_freqs) > 0) {
            p_map[$1"_"$2] = $3;
        }
        close(base_freqs);
    
        print "CHROM", "POS", "LOCAL_SHARING", "EXPECTED_SHARING", "DEVIATION", "VALID_N", "P_BASE", "P_MODERN";
    }
    {
        chr = $1;
        pos = $2;
        p = p_map[chr"_"pos];
    
        if (p == "" || p == 0) next; 
        sum_c = 0;
        sum_c_sq = 0;
        valid_n = 0;
        sum_raw_alleles = 0

        for (i = 3; i <= NF; i++) {

            if ($i ~ /\.\/\./ || $i ~ /\./) continue; 
            split($i, alleles, "/|\\|");
            genotype_score = alleles[1] + alleles[2];
            c = genotype_score - (2 * p);
            sum_c += c;
            sum_c_sq += (c * c);
            valid_n++;
        }

        if (valid_n < 2) next;
        pairs = valid_n * (valid_n - 1);
        local_sharing = (sum_c * sum_c) - sum_c_sq;
        expected_sharing = sum_G_offdiag * (2 * p * (1 - p));
        deviation = (local_sharing - expected_sharing) / pairs;
        p_modern = sum_raw_alleles / (2 * valid_n)
        print chr, pos, local_sharing, expected_sharing, deviation, valid_n, p, p_modern;

    }' $varpos

}

export -f normalise_deltap

# Then stream the variants in vcf format to the function.
# Run in parallel as it is quite slow.
chromosomes=($(bcftools query --format '%CHROM\n' $cultivpedig_ba_bcf | uniq))

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

rm *selscan_af.txt
rm *base_af.txt
rm *varpos.txt

```
## Standard diversity assessments for overlay with previous tests

```bash

# Subset the filtered bcf.
bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --samples-file $pocultivar_ids \
        --force-samples \
        $merge_ba_bcf > $cultivall_ba_bcf

bcftools \
    index \
        --force \
        $cultivall_ba_bcf

bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --samples-file $pedig_all_ids \
        --force-samples \
        $merge_ba_bcf > $pedig_ba_bcf

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

```
# Assess chromosome synteny for swept regions
## Dependencies
minimap

## Code
### Synteny plots for chromosome A08
There are some severe sweeps on chromosome 8.

I will take:
Darmor
zs11
Chikuzen
Chisaya
Haya
Mystic
Trooper
Insignia
Range
Rainbow
X5A

I aim to create a synteny plot to show that all of the cultivars have the same whole chromosome as Haya, Zephyr, Westar and Tower, which are the pedigree founder lines with the X5A haplotype.

```bash

# First, extract all relevant files from the minigraph run.
# If not run already, extract them.
chrom_a8_fastas=(
    $(tar -tf \
        $DATA/workflows/outputs/backup/2026-01-23_14-04-00_cactus.tar.gz  | \
        grep 882 | grep gz | grep -E 'Darmor|zs11|Chikuzen|Chisaya|Haya|Mystic|Trooper|Insignia|Range|Rainbow|X5A')
    )

for fasta in ${chrom_a8_fastas[@]}; do
    tar xvf $DATA/workflows/outputs/backup/2026-01-23_14-04-00_cactus.tar.gz $fasta
done

mv 2026-01-23_14-04-00_cactus/chrom-subproblems/CP123882.1/fasta/*gz $cactus_fastas/a8/ &&
gunzip $cactus_fastas/a8/*gz

# First, run minimap2 for all against X5A.
x5a_fasta=$cactus_fastas/a8/X5A_ref_CP123882.1.fa
mkdir -p $PEDIGREE_OUT_DIR/minimap2/a8
ls $cactus_fastas/a8/*a | grep -v "X5A" | \
        parallel --jobs 13 "minimap2 -cx asm10 $x5a_fasta {} > $PEDIGREE_OUT_DIR/minimap2/a8/{/.}.paf"

```

The dotplots are made in R with pedigree-sel.R.

### Make plots for select varieties with syri
Then, run syri to make an alternative plot. This is to see whether there are any clear<br>
rearrangements in the centromere that are in Zephyr but not Haya / Mystic.
```bash

# Scaffold with ragtag.
conda activate $SOFTWARE/ragtag/ragtag/
mkdir -p $PEDIGREE_OUT_DIR/ragtag

# Make an array of genomes to compare.
GENOMES_COMPARE=(
    $cactus_fastas/a8/Chikuzen_verkko_CP123882.1.fa
    $cactus_fastas/a8/Zephyr_flye_CP123882.1.fa
    $cactus_fastas/a8/Haya_verkko_CP123882.1.fa
    $cactus_fastas/a8/Mystic_verkko_CP123882.1.fa
    $cactus_fastas/a8/Range_verkko_CP123882.1.fa
    $cactus_fastas/a8/Rainbow_flye_CP123882.1.fa
    )

# Scaffold all the genomes with ragtag.
for genome in ${GENOMES_COMPARE[@]}; do

    # Variety name is the file prefix before the first _
    prefix=$(basename $genome | cut -d"_" -f 1)

    # Run ragtag.
    echo "Running ragtag for $prefix"
    ragtag.py \
        scaffold \
        -o $PEDIGREE_OUT_DIR/ragtag/$prefix \
        -u \
        -t 28 \
        $cactus_fastas/a8/X5A_ref_CP123882.1.fa \
        $genome

    # Reformat file.
    awk '$0 ~ />/ {
            if($0 ~ /CP123882/) {sw=1}
                else {sw=0}
            }
         sw == 1 {
            if($0 ~ />/){print ">A08"}
            else{print $0}
            }' $PEDIGREE_OUT_DIR/ragtag/$prefix/ragtag.scaffold.fasta \
            > $PEDIGREE_OUT_DIR/ragtag/${prefix}_scaffold.fasta

done

# Rename the reference chromosomes.
awk '$0 ~ />/ {print ">A08"} $0 !~ />/ {print $0}' \
     $cactus_fastas/a8/X5A_ref_CP123882.1.fa \
     > $cactus_fastas/a8/X5A_ref_CP123882.1_rename.fa

# Run minimap2 then syri for the relevant comparisons.

conda deactivate
conda activate $SOFTWARE/syri/syri
mkdir -p $PEDIGREE_OUT_DIR/syri

GENOMES_A=(
    $PEDIGREE_OUT_DIR/ragtag/Chikuzen_scaffold.fasta
    $PEDIGREE_OUT_DIR/ragtag/Haya_scaffold.fasta
    $PEDIGREE_OUT_DIR/ragtag/Range_scaffold.fasta
    $PEDIGREE_OUT_DIR/ragtag/Mystic_scaffold.fasta
    $PEDIGREE_OUT_DIR/ragtag/Zephyr_scaffold.fasta
    )

GENOMES_B=(
    $PEDIGREE_OUT_DIR/ragtag/Haya_scaffold.fasta
    $PEDIGREE_OUT_DIR/ragtag/Range_scaffold.fasta
    $PEDIGREE_OUT_DIR/ragtag/Mystic_scaffold.fasta
    $PEDIGREE_OUT_DIR/ragtag/Zephyr_scaffold.fasta
    $PEDIGREE_OUT_DIR/ragtag/Rainbow_scaffold.fasta
    )

for ((i=0;i<${#GENOMES_A[@]};i++)); do

    # Prefix based on second genome.
    prefix=$(basename ${GENOMES_B[$i]} | cut -d"_" -f 1)
    echo "Running syri on $prefix"

    # Run minimap2.
    minimap2 \
        -ax asm5 \
        --eqx \
        -t 28 \
        ${GENOMES_A[$i]} \
        ${GENOMES_B[$i]} \
        > $PEDIGREE_OUT_DIR/syri/$prefix.sam
    
    samtools \
        sort \
            --output-fmt SAM \
            $PEDIGREE_OUT_DIR/syri/$prefix.sam \
            > $PEDIGREE_OUT_DIR/syri/$prefix.sorted.sam

    syri \
        -c $PEDIGREE_OUT_DIR/syri/$prefix.sorted.sam \
        -r ${GENOMES_A[$i]} \
        -q ${GENOMES_B[$i]} \
        -k \
        -F S \
        --dir $PEDIGREE_OUT_DIR/syri \
        --prefix ${prefix}_

done

# Generate the config file for plotsr.
printf "%s\n" ${GENOMES_B[@]} | \
    awk 'NR == 1 {print "#file\tname\ttags"}
         {
            split($1,a,"/")
            split(a[length(a)],b,"_")
            print $1"\t"b[1]"\tlw:1.5"
        }' > genomes.txt

syri --nc 5 \
  -c alignment.paf \
  -r X5A_A08.fa \
  -q variety_A08_contig.fa \
  --prefix variety_syri

plotsr \
    --sr $PEDIGREE_OUT_DIR/syri/Haya_syri.out \
    --sr $PEDIGREE_OUT_DIR/syri/Range_syri.out \
    --sr $PEDIGREE_OUT_DIR/syri/Mystic_syri.out \
    --sr $PEDIGREE_OUT_DIR/syri/Zephyr_syri.out \
    --sr $PEDIGREE_OUT_DIR/syri/Rainbow_syri.out \
    --genomes genomes.txt \
    -o $PEDIGREE_OUT_DIR/syri/A08_plotsr.png


```
### Plots for candidate swept regions
There are a few key regions that look like they may have been heavily selected.<br>
The aim here is to pick a couple of key traits, for example erucic acid content,<br>
and show the heavy selection in the candidate region.<br>
#### Fatty acid elongase genes for erucic acid content
I will start with the genes BnaA08g11130D and BnaC03g65980D (Darmor), which are<br>
both involved in erucic acid biosynthesis (being homoeologous).<br>
```bash

# Download the B. napus Darmor annotations.
wget \
    https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/brassica_napus/cds/Brassica_napus.AST_PRJEB5043_v1.cds.all.fa.gz

# Extract them.
zcat Brassica_napus.AST_PRJEB5043_v1.cds.all.fa.gz > $PEDIGREE_IN_DIR/Darmor_annotations.fa

# Extract the known A08 erucic acid variants.
awk '$0 ~ />/ {sw = 0}
     $0 ~ /BnaA08g11130D/ || $0 ~ /BnaC03g65980D/ {sw = 1}
     sw == 1 {print $0}' $PEDIGREE_IN_DIR/Darmor_annotations.fa \
     > $PEDIGREE_IN_DIR/Darmor_EU_genes.fa

# For each genome, run minimap2 to find the best mapping
# sites for the fatty acid elongase genes.
ls $PEDIGREE_IN_DIR/genomes/*gz | \
    parallel --jobs 2 "minimap2 -cx splice:hq --cs {} $PEDIGREE_IN_DIR/Darmor_EU_genes.fa > $PEDIGREE_OUT_DIR/{/.}_EA.paf"

# Run against X5A.
minimap2 \
    -cx splice:hq \
    --cs $ref_fasta \
    $PEDIGREE_IN_DIR/Darmor_EU_genes.fa \
    > $PEDIGREE_OUT_DIR/Darmor_EU_regions.paf

awk '$12 == 60 {print $6"\t"$8"\t"$9}' \
    $PEDIGREE_OUT_DIR/Darmor_EU_regions.paf | \
    sed 's/X5A_ref#0#CP123882.1/A08/g' | \
    sed 's/X5A_ref#0#CP123893.1/C03/g' \
    > $PEDIGREE_OUT_DIR/Darmor_EU_regions.txt

# Extract a 2 kb region surrounding the gene.
for region in $PEDIGREE_OUT_DIR/*_EA.paf; do

    # Get the original fasta name
    fasta="$PEDIGREE_IN_DIR/genomes/$(basename $region | sed 's/_EA.paf//g').gz"

    # Prefix for output files.
    prefix="$(echo $region | sed 's/\.fa.*//g')"

    # Get primary alignments and pad in a BED file.
    padleft=30000
    padright=30000
    awk -v padleft=$padleft \
        -v padright=$padright \
        '$12 == 60 {

            if ($5 == "+") {
                qend = $9 + padright

                if($8 < padleft){
                    qstart = 0

                } else {
                    qstart = $8 - padleft

                }

            } else {
                qstart = $8 - padright; qend = $9 + padleft
                
            }

        printf "%s\t%s\t%s\t%s\t.\t%s\n", $6, qstart, qend, $1, $5
    }' $region > $prefix.bed

    # Extract if necessary.
    if [[ ! -f ${fasta%.*} ]]; then
        zcat $fasta > ${fasta%.*}
    fi

    # Use BED and FASTA to extract the region from the genome.
    bedtools \
        getfasta \
            -fi ${fasta%.*} \
            -fo ${prefix}_EA.fasta \
            -s \
            -name \
            -bed $prefix.bed

    echo "Done for $prefix"

done

# Align A08 regions.
# Get the haplotypes for Haya, Chisaya, Chikuzen, Bronowski, Zephyr, Darmor and 
# cultivars Barossa, Castle, Dunkeld, Eureka, Insignia, Maluka, Marnoo, Mystic, Rainbow,
# Range, Stubby, Shiralee, Trooper, Westar, Yickadee.
VARIETIES_KEEP=(
    Darmorbzh_ref_EA.fasta
    X5A_ref_EA.fasta
    Chisaya_verkko_EA.fasta
    Chikuzen_verkko_EA.fasta
    Haya_verkko_EA.fasta
    Zephyr_flye_EA.fasta
    Westar_flye_EA.fasta     
    Bronowski_flye_EA.fasta
    Drakkar_flye_EA.fasta
    Eureka_verkko_EA.fasta
    Barossa_verkko_EA.fasta
    Yickadee_verkko_EA.fasta
    Maluka_flye_EA.fasta
    Rainbow_flye_EA.fasta
    Mystic_verkko_EA.fasta
    Trooper_verkko_EA.fasta
    Castle_flye_EA.fasta
    Dunkeld_verkko_EA.fasta
    Marnoo_flye_EA.fasta
    Range_verkko_EA.fasta
    Stubby_verkko_EA.fasta
    Insignia_verkko_EA.fasta
)

for each in ${VARIETIES_KEEP[@]}; do
    awk '$1 ~ />/ {sw = 0}
        $1 ~ /CDX76287/ {
            sw = 1; split(FILENAME, a, "/"); split(a[length(a)], b, "_");
            split($1, c, "[>:]"); gname = c[2]; vname = b[1]; print ">"gname";"vname
        }
        $1 !~ />/ && sw == 1 {print $0}' $PEDIGREE_OUT_DIR/$each
done > $PEDIGREE_OUT_DIR/A08_FAE_regions.fasta

awk '$1 ~ />/ {sw = 0}
    $1 ~ /CDX76287/ {sw = 1}
    sw == 1 {print $0}' data/pedigree/inputs/Darmor_EU_genes.fa | \
cat - $PEDIGREE_OUT_DIR/A08_FAE_regions.fasta | \
mafft \
    --auto \
    --thread 20 \
    - \
    > $PEDIGREE_OUT_DIR/A08_FAE_regions_mafft.fasta

# Extract A08 from the biallelic filtered and imputed SNP file.
bcftools \
    view \
        --no-version --output-type b --threads 30 \
        --regions A08 \
        $merge_basnpsmaf_bcf \
        > $a8_basnpsmaf_bcf

bcftools index --force $a8_basnpsmaf_bcf

# Extract A08 from the biallelic filtered and imputed SNP file.
bcftools \
    view \
        --no-version --output-type b --threads 30 \
        --regions C03 \
        $merge_basnpsmaf_bcf \
        > $c3_basnpsmaf_bcf

bcftools index --force $c3_basnpsmaf_bcf

# Get the C07 selection candidates.
minimap2 \
    -cx splice:hq \
    --cs $ref_fasta \
    $PEDIGREE_IN_DIR/c7_selcandidates.fasta \
    > $PEDIGREE_OUT_DIR/c7_selcandidates.paf

awk '$0 ~ /CP123894/ {print $6"\t"$8"\t"$9"\t"$1"\t.\t"$5}' \
    $PEDIGREE_OUT_DIR/c7_selcandidates.paf | \
    sed 's/X5A_ref#0#CP123894.1/C07/g' \
    > $PEDIGREE_OUT_DIR/c7_selcandidates.txt

cut -f 1,2,3 data/pedigree/out/selection_candidates_merged.txt | \
bedtools intersect -a - -b $PEDIGREE_OUT_DIR/c7_selcandidates.txt

```
### Transfer annotations to genomes

## Run pbwt
```bash

# Generate population file.
# Assign each founder to a unique (pseudo) population.
awk '{print $1"\t"$1}' $pbwt_ref_samples > $HAPLOTYPES_IN_DIR/popfile.txt

# Index inputs.
bcftools index --force $pbwt_bcf_ref
bcftools index --force $pbwt_bcf_target_sel

# Get the samples to test.
# This removes samples missing in the final bcf.
bcftools \
    query \
        --list-samples \
        $pbwt_bcf_target_sel > $pbwt_target_samples

# Run pbwt for each chromosome independently.
mkdir -p logs/

chromosomes=($(bcftools query --format '%CHROM\n' $pbwt_bcf_target_sel | uniq))

for chromosome in ${chromosomes[@]}; do

    echo "Preparing input."
    
    # Select ref sites just this chromosome.
    bcftools \
        view \
            --no-version --output-type z --threads 30 \
            --regions $chromosome \
            $pbwt_bcf_ref > $HAPLOTYPES_IN_DIR/pbwtin_ref_$chromosome.vcf.gz

    # Make the sites file.
    bcftools \
        query \
            --format '%POS\n' $HAPLOTYPES_IN_DIR/pbwtin_ref_$chromosome.vcf.gz \
            > $HAPLOTYPES_IN_DIR/sites_$chromosome.txt

    # Make a dummy map.
    echo -e "pd\tgd" > $HAPLOTYPES_IN_DIR/map_$chromosome.txt
    awk '{print $1"\t"1e-8*$1}' \
        $HAPLOTYPES_IN_DIR/sites_$chromosome.txt \
        >> $HAPLOTYPES_IN_DIR/map_$chromosome.txt 

    # Select corresponding sample sites.
    bcftools \
        view \
            --no-version --output-type z --threads 30 \
            --regions $chromosome \
            $pbwt_bcf_target_sel > $HAPLOTYPES_IN_DIR/pbwtin_target_sel_$chromosome.vcf.gz

    # Fix lambda to 10 to represent the average number of meioses based on Cowling 2007.
    lambda=10

    # Fix nmatches to 5 because the reference set is small. If there are more than 5 matches,
    # it will just take the top 5. This is likely better than taking 10 similar matches
    # from the 14 ancestors, since no variety is descended from that many ancestors.
    nmatches=5

    # Lmin is calculated to capture ~1 Mb chunks across the genome. L0 should be a ^2 multiplier.
    lmin=$(awk 'NR > 1 {bp=$1} NR == 2 {start=bp} END {dist=bp-start; printf("%.0f\n", 1000000/(dist/(NR-1)))}' $HAPLOTYPES_IN_DIR/map_$chromosome.txt)
    l0=$((lmin * 4))

    echo "Running for chromosome: $chromosome"
    echo "Lmin: $lmin"
    echo "L0: $l0"
    echo "nmatch: $nmatches"

    # Write reference pbwt files.
    pbwt \
        -readVcfGT $HAPLOTYPES_IN_DIR/pbwtin_ref_$chromosome.vcf.gz \
        -writeAll $HAPLOTYPES_IN_DIR/pbwtin_ref_$chromosome

    # Write samples pbwt files.
    pbwt \
        -readVcfGT $HAPLOTYPES_IN_DIR/pbwtin_target_sel_$chromosome.vcf.gz \
        -writeAll $HAPLOTYPES_IN_DIR/pbwtin_target_sel_$chromosome

    # Find maximum matches between samples and founders.
    pbwt \
        -read $HAPLOTYPES_IN_DIR/pbwtin_ref_$chromosome.pbwt \
        -log logs/pbwt_${chromosome}_log.txt \
        -matchDynamic $HAPLOTYPES_IN_DIR/pbwtin_target_sel_$chromosome.pbwt \
        > $HAPLOTYPES_OUT_DIR/pbwtout_${chromosome}_matchNaive.txt

done

```

# Identify variant near chromosome A08 telomere
Right near the telomere on chromosome A08, there is a big block of non-ref variants only<br>
in cultivars and not in their founder ancestors.<br>
I am susepcting it is a homoeologous reombination, so I am going to look for this.<br>
```bash
# Extract chromosome A08 for the cultivars.
tar xvf \
    $DATA/workflows/outputs/backup/2026-01-23_14-04-00_cactus.tar.gz \
    2026-01-23_14-04-00_cactus/chrom-subproblems/CP123882.1/fasta/Eureka_verkko_CP123882.1.fa.gz \
    2026-01-23_14-04-00_cactus/chrom-subproblems/CP123882.1/fasta/Insignia_verkko_CP123882.1.fa.gz \
    2026-01-23_14-04-00_cactus/chrom-subproblems/CP123882.1/fasta/Range_verkko_CP123882.1.fa.gz \
    2026-01-23_14-04-00_cactus/chrom-subproblems/CP123882.1/fasta/Dunkeld_verkko_CP123882.1.fa.gz


# Then, align each against the full X5A genome.

```
# Identify haplotypes in additional genomes
It will be worth identifying haplotypes in some additional genomes from the sra.<br>
The additional genomes to download are in `doc/info/sra-accessions.txt`.<br>
This is a TSV with SRA accession in column1 and variety ID (to be used in the <br>
bcf) in column2.<br>
```bash

# Path to sra toolkit and kmc.
export PATH="$SOFTWARE/kmc/bin":"$SOFTWARE/sratoolkit/sratoolkit.3.2.1-ubuntu64/bin":$PATH
mkdir -p $PEDIGREE_IN_DIR/sra

# Download the accessions.
while IFS=$'\t' read -r acc varname; do

    echo "Downloading $acc for variety $varname"

    prefetch \
        --ascp-path "/home/ubuntu/.aspera/connect/bin/ascp|/home/ubuntu/.aspera/connect/etc/asperaweb_id_dsa.openssh" \
        $acc

    echo -e "$varname download complete\n"

    echo "Extracting $varname"

    fasterq-dump \
        $acc \
        --outdir $PEDIGREE_IN_DIR/sra \
        --threads 28

    echo -e "$varname extracted\n"

done < doc/info/sra-accessions.txt

awk 'NR > 1 {print $1"\t"$2}' doc/info/chromosomes_rename_key_lens.txt \
    > doc/info/chromosomes_rename_bcftools.txt

# Then run vg giraffe to map the reads.
while IFS=$'\t' read -r acc varname; do

    reads1=$PEDIGREE_IN_DIR/sra/${acc}_1.fastq
    reads2=$PEDIGREE_IN_DIR/sra/${acc}_2.fastq
    mkdir -p ${varname}_kmc_tmp

    echo "Processing $reads1 and $reads2 in ${varname}_kmc_tmp"
    echo -e "${reads1}\\n${reads2}" > ${varname}_kff_input.txt

    kmc \
        -k29 \
        -m20 \
        -t20 \
        -okff \
        @${varname}_kff_input.txt \
        ${varname} \
        ${varname}_kmc_tmp 1> ${varname}_kmc.out 2> ${varname}_kmc.err

    echo -e "K-mer profile for $varname done\n"
    echo "Running vg giraffe for $varname"

    vg giraffe \
        --gbz-name data/graphmap/graph/ccdm-bnpan.gbz \
        --progress \
        --threads 28 \
        --fastq-in $reads1 \
        --fastq-in $reads2 \
        --haplotype-name data/graphmap/graph/ccdm-bnpan.hapl \
        --kff-name ${varname}.kff \
        --index-basename ccdm-bnpan \
        --sample ${varname} 1> ${varname}_giraffe.gam 2> ${varname}_giraffe.log

    echo -e "Finished running vg giraffe for $varname"

    echo "Running vg filter for $varname"

    vg \
        filter \
            --proper-pairs \
            --min-primary 0.9 \
            --frac-score \
            --substitutions \
            --min-end-matches 1 \
            --min-mapq 30 \
            --defray-ends 999 \
            --xg-name $graph_gbz \
            --threads 20 \
            ${varname}_giraffe.gam > ${varname}-filtered.gam

    echo "Running vg pack for $varname"

    vg \
        pack \
            --xg $graph_gbz \
            --gam ${varname}-filtered.gam \
            --packs-out ${varname}.pack \
            --threads 20

    echo -e "Finished running vg pack for $varname"

    echo "Running vg call for $varname"

    vg \
        call \
            --genotype-snarls \
            --all-snarls \
            --sample $varname \
            --pack $varname.pack \
            --ref-sample "X5A_ref" \
            --progress \
            --threads 24 \
            --gbz \
            $graph_gbz > $varname.vcf

    echo -e "Finished running vg call for $varname\n"

    echo "Renaming chromosomes."

    awk '$0 ~ /##contig=/ {split($1,a,"X5A_ref#0#"); print a[1]""a[2]}
    $0 !~ /##contig=/{print $0}' $varname.vcf > ${varname}_rehead.vcf

    bcftools \
        view \
            --threads 30 --no-version --output-type b \
            ${varname}_rehead.vcf \
            > $varname.bcf.gz

    bcftools \
        annotate \
            --threads 30 --no-version --output-type b \
            --rename-chrs doc/info/chromosomes_rename_bcftools.txt \
            $varname.bcf.gz \
            > ${varname}_rename.bcf.gz

    bcftools index --force ${varname}_rename.bcf.gz

    bcftools \
        +setGT ${varname}_rename.bcf.gz \
        -O b -o ${varname}_hetmiss.bcf.gz \
        -- -t q -n "./." --include 'GT == "het"'

    bcftools index ${varname}_hetmiss.bcf.gz

    echo -e "Chromsomes renamed for $varname\n"

done < doc/info/sra-accessions.txt

# https://cdnsciencepub.com/doi/epdf/10.4141/cjps61-028

# I only ran Liho so far. I will merge this with the rest of the genotypes at
# just the called high quality sites.
bcftools \
    query \
        --format '%CHROM\t%POS' \
        $a8_basnpsmaf_bcf \
        > a8_basnpsmaf_pos.txt

bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --regions-file a8_basnpsmaf_pos.txt \
        Liho_hetmiss.bcf.gz > Liho_a8.bcf.gz

bcftools index --force Liho_a8.bcf.gz

bcftools \
    merge \
        --threads 30 --no-version --output-type b \
        $merge_basnpsmaf_bcf Liho_a8.bcf.gz \
        > $PEDIGREE_OUT_DIR/mergeliho_basnpsmaf_a8.bcf.gz

bcftools \
    query \
        --format '%CHROM\t%POS' \
        $c3_basnpsmaf_bcf \
        > c3_basnpsmaf_pos.txt

bcftools \
    view \
        --threads 30 --no-version --output-type b \
        --regions-file c3_basnpsmaf_pos.txt \
        Liho_hetmiss.bcf.gz > Liho_c3.bcf.gz

bcftools index --force Liho_c3.bcf.gz

bcftools \
    merge \
        --threads 30 --no-version --output-type b \
        $c3_basnpsmaf_bcf Liho_c3.bcf.gz \
        > $PEDIGREE_OUT_DIR/mergeliho_basnpsmaf_c3.bcf.gz
```