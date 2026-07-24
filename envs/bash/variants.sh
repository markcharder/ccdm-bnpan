#!/usr/bin/env bash

WF_DIR=/usr/local/share/mount_b/202411_brassica_napus_pangenome_metaqtls/mark/workflows/ccdm-bnpan

# Set $SOFTWARE to local siftware install directory.
export PATH=$SOFTWARE/plink/plink2:$SOFTWARE/bcftools/plugins/bcftools/:$SOFTWARE/gcta/gcta-1.95.1-linux-x86_64:$PATH
export BCFTOOLS_PLUGINS=$SOFTWARE/bcftools/plugins/bcftools/plugins/

# Input directories.
VARIANTS_DIR="$WF_DIR/data/variants"
GRAPH_DIR="$VARIANTS_DIR/graph/"
READS_DIR="$VARIANTS_DIR/reads/"
MERGE_DIR="$VARIANTS_DIR/merged/"
SUBSET_DIR="$VARIANTS_DIR/subset/"
KINSHIP_DIR="$VARIANTS_DIR/kinship/"
REF_DIR="$VARIANTS_DIR/reference/"

# Set reference name and FASTA.
ref_name=X5A_ref
ref_fasta=$REF_DIR/X5A_ref.fa

# Capture reference index files.
ref_idx=$ref_fasta.fai
ref_idx_simple=$ref_fasta.simple.fai
ref_idx_chromids=$ref_fasta.chrids.fai

# A file for renaming chromosomes.
chromids_file=$REF_DIR/chromids.txt
chromids_pansn_file=$REF_DIR/chromids_pansn.txt

ncbiaccs=($(cat $ref_idx | cut -f 1 | sed 's/X5A_ref#0#//g'))
chromids=($(echo A0{1..9} A10 C0{1..9}))

# Make files with and without the panSN prefix.
for ((i=0;i<${#ncbiaccs[@]};i++)); do
    echo -e "${ncbiaccs[$i]}\t${chromids[$i]}"
done > $chromids_file

for ((i=0;i<${#ncbiaccs[@]};i++)); do
    echo -e "X5A_ref#0#${ncbiaccs[$i]}\t${chromids[$i]}"
done > $chromids_pansn_file

# Make an fai to update the header.
sed 's/X5A_ref#0#//' $ref_idx > $ref_idx_simple
cut -f 2 $chromids_file | \
    paste - $ref_idx_simple | \
    cut -f 1,3,4,5,6 > $ref_idx_chromids

graph_gbz=$GRAPH_DIR/ccdm-bnpan.gbz
graph_vcf=$GRAPH_DIR/ccdm-bnpan.vcf
graph_bcf=$GRAPH_DIR/ccdm-bnpan.bcf.gz
graph_snarls=$GRAPH_DIR/ccdm-bnpan.snarls

reads_vcf="$READS_DIR/bcftools-merged_f.vcf.gz"
reads_filt_bcf="$READS_DIR/bcftools-merged_f.filtered.bcf.gz"

merge_full_bcf="$MERGE_DIR/merge_full.bcf.gz"
merge_filt_bcf="$MERGE_DIR/merge_filt.bcf.gz"
merge_filtsub_bcf="$MERGE_DIR/merge_filtsub.bcf.gz"
merge_ba_bcf="$MERGE_DIR/merge_ba.bcf.gz"
merge_basnps_bcf="$MERGE_DIR/merge_basnps.bcf.gz"
merge_basvs_bcf="$VARIANTS_DIR/merged/merge_basvs.bcf.gz"

merge_basnpsmaf_bcf="$MERGE_DIR/merge_basnpsmaf.bcf.gz"
merge_basnpsmaf_vcf="$SUBSET_DIR/merge_basnpsmaf.vcf.gz"
merge_basnpsmaf_impute_vcf="$SUBSET_DIR/merge_basnpsmaf_impute"
merge_basnpsmaf_impute_bcf="$SUBSET_DIR/merge_basnpsmaf_impute.bcf.gz"

merge_ldfiltpca_bcf="$MERGE_DIR/merge_ldfiltpca.bcf.gz"
merge_ldfiltpca_plink="$KINSHIP_DIR/merge_ldfiltpca_plink"
merge_ldfiltpca_grm="$KINSHIP_DIR/merge_ldfiltpca_grm"
merge_ldfiltgp_bcf="$MERGE_DIR/merge_ldfiltgp.bcf.gz"
merge_ldfiltgp_plink="$KINSHIP_DIR/merge_ldfiltgp_plink"
merge_ldfiltgp_grm="$KINSHIP_DIR/merge_ldfiltgp_grm"

merge_basnpsmaf_impute_sp_bcf="$SUBSET_DIR/merge_basnpsmaf_impute_sp.bcf.gz"
merge_basnpsmaf_impute_wi_bcf="$SUBSET_DIR/merge_basnpsmaf_impute_wi.bcf.gz"
merge_ldfiltpca_sp_bcf="$SUBSET_DIR/merge_ldfiltpca_sp.bcf.gz"
merge_ldfiltpca_wi_bcf="$SUBSET_DIR/merge_ldfiltpca_wi.bcf.gz"
merge_ldfiltgp_sp_bcf="$SUBSET_DIR/merge_ldfiltgp_sp.bcf.gz"
merge_ldfiltgp_wi_bcf="$SUBSET_DIR/merge_ldfiltgp_wi.bcf.gz"
merge_ldfiltgp_sp_kin="$KINSHIP_DIR/merge_ldfiltgp_sp_kin"
merge_sp_ldfiltgp_plink="$KINSHIP_DIR/merge_sp_ldfiltgp_plink"
merge_wi_ldfiltgp_plink="$KINSHIP_DIR/merge_wi_ldfiltgp_plink"

merge_basnpsmaf_plink="$KINSHIP_DIR/merge_basnpsmaf_plink"
merge_filtsub_spring_bcf="$SUBSET_DIR/merge_filtsub_spring.bcf.gz"
merge_filtsub_winter_bcf="$SUBSET_DIR/merge_filtsub_winter.bcf.gz"
merge_ba_spring_bcf="$SUBSET_DIR/merge_ba_spring.bcf.gz"
merge_ba_winter_bcf="$SUBSET_DIR/merge_ba_winter.bcf.gz"
merge_basnps_spring_bcf="$SUBSET_DIR/merge_basnps_spring.bcf.gz"
merge_basnps_winter_bcf="$SUBSET_DIR/merge_basnps_winter.bcf.gz"
merge_basnpsmaf_spring_bcf="$SUBSET_DIR/merge_basnpsmaf_spring.bcf.gz"
merge_basnpsmaf_winter_bcf="$SUBSET_DIR/merge_basnpsmaf_winter.bcf.gz"
merge_ldfilt_spring_bcf="$SUBSET_DIR/merge_ldfilt_spring.bcf.gz"
merge_ldfilt_winter_bcf="$SUBSET_DIR/merge_ldfilt_winter.bcf.gz"

pedig_basnpsmaf_bcf="$SUBSET_DIR/pedig_basnpsmaf.bcf.gz"
pedig_ldfilt_bcf="$SUBSET_DIR/pedig_ldfilt.bcf.gz"
pedig_ldfilt_nz="$SUBSET_DIR/pedig_ldfilt_nz.bcf.gz"
pedig_basnpsmaf_plink="$KINSHIP_DIR/pedig_basnpsmaf_plink"
pedig_ldfiltpca_plink="$KINSHIP_DIR/pedig_ldfiltpca_plink"
pedig_base_af="$KINSHIP_DIR/pedig_base_af.txt"

cultiv_basnpsmaf_bcf="$SUBSET_DIR/cultiv_basnpsmaf.bcf.gz"
cultiv_basnpsmaf_plink="$KINSHIP_DIR/cultiv_basnpsmaf_plink"
cultiv_ldfilt_bcf="$SUBSET_DIR/cultiv_ldfilt.bcf.gz"
cultiv_ldfilt_kin="$KINSHIP_DIR/cultiv_ldfilt_kin"

merge_basnpsmaf_spring_kin="$KINSHIP_DIR/merge_basnpsmaf_spring_kin"
merge_basnpsmaf_winter_kin="$KINSHIP_DIR/merge_basnpsmaf_winter_kin"
merge_ldfilt_spring_kin="$KINSHIP_DIR/merge_ldfilt_spring_kin"
merge_ldfilt_winter_kin="$KINSHIP_DIR/merge_ldfilt_winter_kin"

founders_af="$SUBSET_DIR/founders_af.txt"
founders_winter_af="$SUBSET_DIR/founders_winter_af.txt"
founders_spring_af="$SUBSET_DIR/founders_spring_af.txt"
founders_semiwinter_af="$SUBSET_DIR/founders_semiwinter_af.txt"
founders_zero="$SUBSET_DIR/founders_zero.txt"
founders_nz="$SUBSET_DIR/founders_nz.txt"
founders_all_af="$SUBSET_DIR/founders_all_af.txt"
