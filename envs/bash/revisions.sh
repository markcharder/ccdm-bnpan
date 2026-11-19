
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

# Set $SOFTWARE to local siftware install directory.
export PATH="$SOFTWARE/plink/plink2":"$SOFTWARE/bcftools/plugins/bcftools/":"$SOFTWARE/gcta/gcta-1.95.1-linux-x86_64":"$PATH"
export BCFTOOLS_PLUGINS="$SOFTWARE/bcftools/plugins/bcftools/plugins/"
export PATH=$PATH:"$SOFTWARE/selscan/selscan/bin/linux/"

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

cultivpedig_ba_bcf="$PEDIGREE_IN_DIR/cultivpedig_ba.bcf.gz"
founders_ba_bcf="$PEDIGREE_IN_DIR/founders_ba.bcf.gz"
founders_winter_af_full="$PEDIGREE_IN_DIR/founders_winter_af_full.txt"
founders_semiwinter_af_full="$PEDIGREE_IN_DIR/founders_semiwinter_af_full.txt"
founders_spring_af_full="$PEDIGREE_IN_DIR/founders_spring_af_full.txt"
founders_all_af_full="$PEDIGREE_IN_DIR/founders_all_af_full.txt"
pedig_base_af_full="$PEDIGREE_IN_DIR/pedig_base_af_full.txt"
selscan_deltap="$PEDIGREE_IN_DIR/selscan_deltap"

pocultivar_ids=$PEDIGREE_IN_DIR/pocultivar_ids.txt
pocultivar_basnpsmaf_impute_vcf=$PEDIGREE_IN_DIR/pocultivar_ldfiltsel.vcf.gz
pocultivar_basnpsmaf_impute_bcf=$PEDIGREE_IN_DIR/pocultivar_basnpsmaf_impute.bcf.gz
pedig_all_basnpsmaf_impute_bcf=$PEDIGREE_IN_DIR/pedig_all_basnpsmaf_impute.bcf.gz

# This is added - run the iHS scan and other diversity tests also on
# Australian breeding lines.
ausbl_ids=$PEDIGREE_IN_DIR/ausbl_ids.txt
ausbl_basnpsmaf_impute_vcf=$PEDIGREE_IN_DIR/ausbl_basnpsmaf_impute.vcf.gz
ausbl_basnpsmaf_impute_bcf=$PEDIGREE_IN_DIR/ausbl_basnpsmaf_impute.bcf.gz
cultivall_ba_bcf=$PEDIGREE_IN_DIR/cultivall_ba.bcf.gz
pedig_ba_bcf=$PEDIGREE_IN_DIR/pedig_ba.bcf.gz
cultivall_ba_pi="$PEDIGREE_OUT_DIR/cultivall_ba_pi.txt"
cultivall_ba_tajima="$PEDIGREE_IN_DIR/cultivall_ba_tajima"
acultivar_ids=$PEDIGREE_IN_DIR/acultivar_ids.txt
acultivar_basnpsmaf_bcf=$PEDIGREE_IN_DIR/acultivar_basnpsmaf.bcf.gz
