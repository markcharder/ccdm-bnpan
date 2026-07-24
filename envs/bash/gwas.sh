source "$(dirname ${BASH_SOURCE[0]})/variants.sh"

# GWAS IO directories.
GWAS_IN_DIR=$WF_DIR/data/gwas/inputs/
GWAS_OUT_DIR=$WF_DIR/data/gwas/out/

# Pre-existing input files other than those in variants.sh.
fta_ids=$GWAS_IN_DIR/fta_ids.txt
glasshouse_ids=$GWAS_IN_DIR/glasshouse_ids.txt
fta_nw_ids=$GWAS_IN_DIR/fta_nw_ids.txt

# Prepared variants for GWAS.
merge_basnpsmaf_iprune_bcf=$GWAS_IN_DIR/merge_basnpsmaf_iprune.bcf.gz
merge_basnpsmaf_iprune_plink=$GWAS_IN_DIR/merge_basnpsmaf_iprune_plink
merge_basnpsmaf_iprune_snpmat=$GWAS_IN_DIR/merge_basnpsmaf_iprune.dat
merge_basnpsmaf_iprune_map=$GWAS_IN_DIR/merge_basnpsmaf_iprune.map
merge_basnpsmaf_iprune_ids=$GWAS_IN_DIR/merge_basnpsmaf_iprune.id

# Prepare relationship matrix for PCA.
merge_basnpsmaf_ldfilt_bcf=$GWAS_IN_DIR/merge_basnpsmaf_ldfilt.bcf.gz
merge_basnpsmaf_ldfilt_plink=$GWAS_IN_DIR/merge_basnpsmaf_ldfilt_plink
merge_basnpsmaf_ldfilt_grm=$GWAS_IN_DIR/merge_basnpsmaf_ldfilt_grm
