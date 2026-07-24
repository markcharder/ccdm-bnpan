# Hacks
## Renaming chromosomes
I have discovered just now that all the chromosomes I have labelled are wrong. I got the X5A chromosome IDs mixed up. I have created a file this way for updating chromosome names.
```bash

# First, extract all chromosome fastas from the minigraph cactus alignment.
ref_chrom_fastas=(
    $(tar -tf $DATA/workflows/outputs/backup/2026-01-23_14-04-00_cactus.tar.gz | \
        grep fasta | grep gz | grep -v MINIGRAPH | grep _ref)
    )

mkdir -p tmp && cd tmp

for fasta in ${ref_chrom_fastas[@]}; do
    tar xvf $DATA/workflows/outputs/backup/2026-01-23_14-04-00_cactus.tar.gz $fasta
done

for fasta in ${ref_chrom_fastas[@]}; do
    if [[ $fasta != *AMBIGUOUS* ]]; then
        suffix=${fasta##*/}
        acc=$(echo ${suffix%.fa*} | cut -d"_" -f 3)
        zcat $fasta | grep ">" | \
            awk -F'[|_]' -v acc=$acc '$0 !~ /scaffold/ && $0 !~ /X5A/ {
                print acc"\t"$NF
            }'
    fi
done | uniq | sort -V -k1,1 > real_chrom_ids.txt

sort -V -k1,1 ../data/variants/reference/chromids.txt > fake_chrom_ids.txt
paste real_chrom_ids.txt fake_chrom_ids.txt | \
    cut -f 1,2,4 | awk 'BEGIN {print "ACCESSION\tREAL_ID\tFAKE_ID"} {print $0}' \
    > ../doc/info/chromosomes_rename_key.txt

```
# Pan-genome graph construction with cactus
## Dependencies
cactus v3.1.2
## References
X5A genome https://doi.org/10.1093/hr/uhad171\n
Cactus algorithm https://doi.org/10.1038/s41586-020-2871-y
## Background
Using the recent high quality X5A genome as the reference, minigraph cactus will be run on a few high quality Brassica napus genomes and the new genomes generated in our study. This assumes all necessary assemblies are in 'data/'. Information on the assemblies is in 'info/assemblies.csv'.\n
## Code
The input for cactus-pangenome is a config file with variety name and genome FASTA file. Assuming we are 1) in the root of the github repository and all genomes are in data/ and 2) the files are named as '<sample>.fasta.gz' and gzipped (would also work without). We do the following:
```bash
mkdir -p conf
for file in `pwd`/data/*gz; do
    sample=$(echo $(basename $file) | cut -d"." -f 1)
    echo -e "$sample\t$file"
done > conf/cactus_genomes.txt
```
Then run cactus. In this case, cactus is installed in ~/.local/software/cactus-bin-v3.1.2 based on the instructions on the cactus github page (https://github.com/ComparativeGenomicsToolkit/cactus.git). Note, the --refContigs option should have contig names corresponding to the assembled chromosomes in the chosen reference, which in our case is X5A_ref (in the file X5A_ref.fasta.gz).
```bash
mkdir -p .cactus_tmp/work
mkdir -p results/
source ~/.local/software/cactus-bin-v3.1.2/venv-cactus-v3.1.2/bin/activate
TIMESTAMP="$(date +"%Y-%m-%d_%H-%M-%S")"
cactus-pangenome \
    --refContigs A0{1..9} A10 C0{1..9} \
    --outDir results/${TIMESTAMP}_cactus \
    --outName ccdm-bnpan \
    --reference X5A_ref \
    --maxCores 30 \
    --maxMemory 200Gi \
    --gbz clip \
    --gfa clip \
    --haplo clip \
    --logDebug \
    --logFile .cactus_tmp/${TIMESTAMP}_cactus.log \
    --workDir `pwd`/.cactus_tmp/work/  \
    --cleanWorkDir onSuccess \
    --stopOnFirstFailure 1 \
    .jobStore-${TIMESTAMP}/ \
    conf/cactus_genomes.txt 1> .cactus_tmp/${TIMESTAMP}_cactus.stdout 2> .cactus_tmp/${TIMESTAMP}_cactus.stderr

```
W       Maluka_flye     0       contig_6818     167     45513
# Identification of telomeres
## Dependencies
tidk 0.2.7
## References
GNU parallel https://doi.org/10.5281/zenodo.17692695
tidk software https://doi.org/10.1093/bioinformatics/btaf049
Telomere motif used matches this paper https://doi.org/10.1093/hr/uhad171
## Background
To determine how contiguous assemblies are, primary contigs will be assessed for the presence of telomers at either end. The telomere motif used is 'TTTAGGG'.
## Code
Assuming the primary assemblies are in ccdm-bnpan/data, run the following code to search for telomeres.
```bash
mkdir -p results/tidk/ logs/
WD=`pwd`
ls `pwd`/data/assemblies/* | \
    parallel --jobs 10 "tidk search --string TTTAGGG --output {/.} --dir $WD/results/tidk/ {} 1> $WD/logs/tidk_{/.}.log 2>&1"
```
# Downloading and clustering fastqs
## Dependencies
sratoolkit 3.2.1 (prefetch)
sratoolkit 3.3.0 (fasterq-dump)
## References
## Background
There are a lot of Illumina reads files in the NCBI sequence read archive for Brassica napus. Many of them are for the same varieties that we sequenced. I aim to download a significant number of these reads files. To reduce redundancy, I was thinking of retaining the file with the most coverage for any almost identical varieties. For varieties that we sequenced that already have data, I will use coverage from both files if they are almost identical.
## Code
Download accessions first from the SRA. This assumes sratoolkit is in $PATH.
```bash
ASCP=/home/ubuntu/.aspera/

cut -d, -f 9 doc/info/sample-ids.csv | grep -vE "NONE|SEQUENCING_ACCESSION" > doc/info/accessions.txt
grep "^SRR" doc/info/accessions.txt | sort -u > doc/info/accessions-sra.txt
mkdir -p data/illumina/

cd data/illumina/

cat ../../doc/info/accessions-sra.txt | \
    parallel --jobs 4 "prefetch --resume yes --ascp-path '$ASCP/connect/bin/ascp|$ASCP/connect/etc/asperaweb_id_dsa.openssh' --progress {}"

# Extract with fasterq-dump.
ls -d SRR* | parallel --jobs 10 "fasterq-dump --split-files {}"
```
Then, download files from the Chinese genome sequence archive.
```bash
# Run this first. There are 35 accessions.
# Find their download paths at https://ngdc.cncb.ac.cn/gsa/browse/
grep -v "^SRR" doc/info/accessions.txt  | sort -u > doc/info/accessions-gsa.txt

# Once download paths have been added to this file, replacing accessions, run this.
cd data/illumina/
cat ../../doc/info/accessions-gsa.txt | parallel --jobs 10 "wget {}"
```
# Renaming contigs in reference genomes to match annotation outputs
```bash
mkdir -p data/raw/genetests/
ref_assemblies=($(ls data/raw/assemprep/assemblies/*gz | grep -v "X5A"))
printf "%s\n" ${ref_assemblies[@]} | parallel --jobs 8 "zcat {} > tmp/{/.}"
ref_assemblies_ex=($(ls tmp/*.fa))

for assem in ${ref_assemblies_ex[@]}; do

    suffix=${assem##*/}
    prefix=${suffix%.*}
    scaf_prefix=${prefix}_0_

    if [ $suffix == 'bnapus_Darmorbzh.fa' ]; then

        awk -v sp=$scaf_prefix '$0 ~ />/ {
            if($0 !~ /scaffold/) {print ">"sp""substr($1,2)}
            else {print $0}
        } $0 !~ />/ {print $0}' $assem > tmp/${prefix}_rn.fasta

    elif [ $suffix == 'bnapus_zs11.fa' ]; then

        awk -v sp=$scaf_prefix '$0 ~ />/ {
            if($0 ~ /[AC](0[0-9]|10)/) {split($0,a,"scaffold"); print ">"sp""a[2]}
            else {print $0}
        } $0 !~ />/ {print $0}' $assem > tmp/${prefix}_rn.fasta

    else

        awk -v sp=$scaf_prefix '$0 ~ />/ {
        if($0 !~ /scaffold/) {split($0,a,"_"); print ">"sp""a[2]}
        else {print $0}
    } $0 !~ />/ {print $0}' $assem > tmp/${prefix}_rn.fasta

    fi

done

```

# SANDPIT
# Call variants from the graph without short reads
## Dependencies
vg 1.70.0
bgzip 1.11
bcftools 1.19
## Code
Using vg deconstruct, get the variants from the graph.
```bash

export PATH=$PATH:$SOFTWARE/plink/plink2
ref_name=X5A_ref
ref_dir=data/graphcall/reference/
ref_fasta=$ref_dir/X5A_ref.fa
ref_idx=$ref_fasta.fai
ref_idx_simple=$ref_fasta.simple.fai
ref_idx_chromids=$ref_fasta.chrids.fai
chromids_file=$ref_dir/chromids.txt
chromids_pansn_file=$ref_dir/chromids_pansn.txt
graph_dir=data/deconstruct/graph/
graph_gbz=$graph_dir/ccdm-bnpan.gbz
graph_vcf=$graph_dir/ccdm-bnpan.vcf
graph_bcf=$graph_dir/ccdm-bnpan.bcf.gz
graph_snarls=$graph_dir/ccdm-bnpan.snarls
gwas_dir=data/gwas
gwas_inputs=$gwas_dir/inputs/
gwas_out=$gwas_dir/out
gwas_bcfmerge=$gwas_inputs/graph_reads_merge.bcf.gz
gwas_bcf_full=$gwas_inputs/gwasin_full.bcf.gz
gwas_bcf_snps=$gwas_inputs/gwasin_snps.bcf.gz
gwas_bcf_svs=$gwas_inputs/gwasin_svs.bcf.gz
gwas_bcf_nsvs=$gwas_inputs/gwasin_nsvs.bcf.gz
gwas_vcf_snps=$gwas_inputs/gwasin_snps.vcf.gz
gwas_vcf_full=$gwas_inputs/gwasin_full.vcf.gz
gwas_vcf_impute=$gwas_inputs/gwasin_impute
gwas_bcf_ldfilt=$gwas_inputs/gwasin_ldfilt
gwas_bcf_impute=$gwas_inputs/gwasin_impute
gwas_kinship=$gwas_inputs/gwasin_kinship
gwas_pca=$gwas_inputs/gwasin_pca
gwas_phenos=$gwas_inputs/phenos_eigenvec.txt
graphcall_dir=data/graphcall/out/bcftools_merge/
graphcall_merge=$graphcall_dir/bcftools-merged_f.vcf.gz
graphcall_merge_filt=${graphcall_merge%.vcf*}.filtered.bcf.gz

export PATH="$SOFTWARE"/bcftools:$PATH
export BCFTOOLS_PLUGINS="$SOFTWARE"/bcftools/plugins/bcftools/plugins

mkdir -p $gwas_inputs
mkdir -p $gwas_out

# A file for renaming chromosomes.
ncbiaccs=($(cat $ref_idx | cut -f 1 | sed 's/X5A_ref#0#//g'))
chromids=($(echo A0{1..9} A10 C0{1..9}))

# Make files with and without the panSN prefix.
for ((i=0;i<${#ncbiaccs[@]};i++)); do
    echo -e "${ncbiaccs[$i]}\t${chromids[$i]}"
done > data/graphcall/reference/chromids.txt

for ((i=0;i<${#ncbiaccs[@]};i++)); do
    echo -e "X5A_ref#0#${ncbiaccs[$i]}\t${chromids[$i]}"
done > data/graphcall/reference/chromids_pansn.txt

# Make an fai to update the header.
sed 's/X5A_ref#0#//' $ref_idx > $ref_idx_simple
cut -f 2 $chromids_file | \
    paste - $ref_idx_simple | \
    cut -f 1,3,4,5,6 > $ref_idx_chromids

# Call variants from the graph only with deconstruct. #
vg \
    deconstruct \
        --all-snarls \
        --path-prefix $ref_name \
        --snarls $graph_snarls \
        --threads 30 \
        $graph_gbz > $graph_vcf

# Set conflicting genotypes to missing, rename #
# chromosomes, remove INFO except AT and all FILTER IDs. #
./util/reannotate_vcf_conflicts.py \
    $graph_vcf | \
bcftools \
    annotate \
        --no-version \
        --include 'LV=0' \
        --remove INFO/AC,INFO/AF,INFO/NS,INFO/AN,FILTER,INFO/LV,INFO/PS,FILTER \
        --output-type b \
        --output $graph_bcf \
        --threads 30 \
        --rename-chrs $chromids_pansn_file \
        -

# Set heterozygous, low quality (GQ < 20 + DP < 2) genotypes to missing.
# Rename chromosomes. #
bcftools \
    +setGT \
        --threads 30 \
        --output-type u \
        $graphcall_merge \
        -- \
        --include 'GT == "het" | FMT/GQ < 20 | FMT/DP < 2' \
        --target-gt q \
        --new-gt './.' | \
bcftools \
    annotate \
        --no-version \
        --remove INFO/LV,INFO/PS,INFO/DP,FILTER,FORMAT/AD,FORMAT/MAD,FORMAT/DP,FORMAT/GL,FORMAT/GQ,FORMAT/GP,FORMAT/XD \
        --output-type b \
        --output $graphcall_merge_filt \
        --threads 30 \
        --rename-chrs $chromids_file \
        -

# Index the prepared VCFs. #
bcftools \
    index \
    --threads 30 \
    $graph_bcf

bcftools \
    index \
    --threads 30 \
    $graphcall_merge_filt

# Merge prepared VCFs. #
bcftools \
    merge \
    --threads 30 \
    --output-type b \
    --output $gwas_bcfmerge \
    $graph_bcf $graphcall_merge_filt

# Make a file for standard GWAS containing all types of #
# biallelic variants, filtered for MAF (minor allele frequency) #
# and fraction missing (below 10%). #
bcftools \
    filter \
        --threads 30 \
        --include 'F_MISSING < 0.2' \
        --output-type b \
        $gwas_bcfmerge > $gwas_bcf_full
bcftools \
    filter \
        --threads 30 \
        --include 'TYPE = "snp"' \
        --output-type b \
        $gwas_bcf_full > $gwas_bcf_snps
bcftools \
    filter \
    --threads 30 \
    --include 'ILEN > 49 || ILEN < -49' \
    --output-type b \
    $gwas_bcf_full > $gwas_bcf_svs
bcftools \
    filter \
    --threads 30 \
    --exclude 'ILEN > 49 || ILEN < -49' \
    --output-type b \
    $gwas_bcf_full > $gwas_bcf_nsvs

# LD prune for PCA #
bcftools \
    filter \
    --include 'MAF >= 0.05 && N_ALT == 1' \
    --output-type u \
    --threads 30 \
    $gwas_bcf_snps | \
bcftools \
    +prune \
        --max r2=0.6 \
        --output-type b \
        - > $gwas_bcf_ldfilt.snps.bcf.gz

# Make kinship matrices. #
plink2 \
    --bcf $gwas_bcf_ldfilt.snps.bcf.gz \
    --make-rel square \
    --out $gwas_kinship.ldfiltsnps
plink2 \
    --bcf $gwas_bcf_snps \
    --make-rel square \
    --out $gwas_kinship.snps
plink2 \
    --bcf $gwas_bcf_nsvs \
    --make-rel square \
    --out $gwas_kinship.nsv
plink2 \
    --bcf $gwas_bcf_svs \
    --make-rel square \
    --out $gwas_kinship.sv

# Prepare input for BLINK. #
# Impute missing alleles with beagle. #
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


# Diagnosing kmer distributions for reads
## Dependencies
seqkit 2.3.0
kmc 3.2.4
## Code
I have two illumina sequence files for testing in 'kmer_test/'. They are forward reads from our study '10_227CTHLT4_GACTGAGTAG-GTTGATAGTG_L001_R1_trimmed.fastq.gz' and reads from a previous study 'CRR045738_f1.fq.gz'.

File sizes are:
9.4G 10_227CTHLT4_GACTGAGTAG-GTTGATAGTG_L001_R1_trimmed.fastq.gz
1.6G CRR045738_f1.fq.gz

First, I will estimate read coverage, based on a Brassica napus genome size of 1.2 Gb.
```bash
cd kmer_test

# Get read stats.
seqkit \
    stats \
        --all \
        10_227CTHLT4_GACTGAGTAG-GTTGATAGTG_L001_R1_trimmed.fastq.gz > ccdm-reads-stats.txt

seqkit \
    stats \
        --all \
        CRR045738_f1.fq.gz > published-reads-stats.txt

# Make kmer histograms.
mkdir -p .kmc_tmp/

kmc \
    -k29 \
    -m40 \
    -t20 \
    CRR045738_f1.fq.gz \
    published-reads-k29 \
    .kmc_tmp/

kmc \
    -k31 \
    -m40 \
    -t20 \
    CRR045738_f1.fq.gz \
    published-reads-k31 \
    .kmc_tmp/

kmc \
    -k29 \
    -m40 \
    -t20 \
    10_227CTHLT4_GACTGAGTAG-GTTGATAGTG_L001_R1_trimmed.fastq.gz \
    ccdm-reads-k29 \
    .kmc_tmp/

kmc \
    -k31 \
    -m40 \
    -t20 \
    10_227CTHLT4_GACTGAGTAG-GTTGATAGTG_L001_R1_trimmed.fastq.gz \
    ccdm-reads-k31 \
    .kmc_tmp/

# Make histograms.
kmc_tools \
    transform \
    published-reads-k29 \
    histogram \
    published-reads-k29.hist

kmc_tools \
    transform \
    published-reads-k31 \
    histogram \
    published-reads-k31.hist

kmc_tools \
    transform \
    ccdm-reads-k29 \
    histogram \
    ccdm-reads-k29.hist

kmc_tools \
    transform \
    ccdm-reads-k31 \
    histogram \
    ccdm-reads-k31.hist

# Make plots.
for file in *.hist; do
    python \
        ../bin/plot_kmc_histogram.py \
        $file -o ${file%.*}.png
done

```
# Clustering fastqs
The first step is to cluster the downloaded FASTQ files to get groups of highly similar samples.
```bash
mash sketch -r -m 2 -k 21 -s 10000 -o all_samples -l file_list.txt
```
# Testing vcfbub
```bash
ls *vcf | parallel --jobs 2 "wc -l {} && sed 's/X5A_ref#0#//g' {} > {.}-rf.vcf"
ls *-rf.vcf | parallel --jobs 2 "time vcfbub --max-level 0 --input {} > {.}-lv0.vcf"
ls *-rf.vcf | parallel --jobs 2 "wc -l {}"
ls *-lv0.vcf | parallel --jobs 2 "wc -l {}"
for vcf in *-lv0.vcf; do
    bgzip --threads 32 $vcf
done
for vcfz in *gz; do
    bcftools index --threads 32 $vcfz
done
ls *vcf.gz > file-list.txt
bcftools merge --file-list file-list.txt --threads 32 --output-type z --output test-merge.vcf.gz
```
