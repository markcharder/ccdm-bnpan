source "$(dirname ${BASH_SOURCE[0]})/variants.sh"

# Set $SOFTWARE to local software install directory.
export PATH=$SOFTWARE/selscan/selscan/bin/linux/:$SOFTWARE/bcftools/plugins/bcftools/:$SOFTWARE/plink/plink2:$SOFTWARE/pbwt/SparsePainter/:$SOFTWARE/pbwt/pbwt/:$PATH
export BCFTOOLS_PLUGINS=$SOFTWARE/bcftools/plugins/bcftools/plugins/
WF_DIR=/usr/local/share/mount_b/202411_brassica_napus_pangenome_metaqtls/mark/workflows/ccdm-bnpan

# Input / output directories.
PEDIGREE_IN_DIR=$WF_DIR/data/pedigree/inputs/
PEDIGREE_OUT_DIR=$WF_DIR/data/pedigree/out/

# Estimating diversity loss from founders to cultivars
pedig_founder_ids=$PEDIGREE_IN_DIR/pedig_founder_ids.txt
pedig_cultivar_ids=$PEDIGREE_IN_DIR/pedig_cultivar_ids.txt
other_cultivar_ids=$PEDIGREE_IN_DIR/other_cultivar_ids.txt
pocultivar_ids=$PEDIGREE_IN_DIR/pocultivar_ids.txt
pfounder_basnpsmaf_bcf=$PEDIGREE_IN_DIR/pfounder_basnpsmaf.bcf.gz
pcultivar_basnpsmaf_bcf=$PEDIGREE_IN_DIR/pcultivar_basnpsmaf.bcf.gz
ocultivar_basnpsmaf_bcf=$PEDIGREE_IN_DIR/ocultivar_basnpsmaf.bcf.gz
chromosomes=(A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 C01 C02 C03 C04 C05 C06 C07 C08 C09)
pfounder_basvs_bcf=$PEDIGREE_IN_DIR/pfounder_basvs.bcf.gz
pcultivar_basvs_bcf=$PEDIGREE_IN_DIR/pcultivar_basvs.bcf.gz
ocultivar_basvs_bcf=$PEDIGREE_IN_DIR/ocultivar_basvs.bcf.gz
acultivar_ids=$PEDIGREE_IN_DIR/acultivar_ids.txt
acultivar_basnpsmaf_bcf=$PEDIGREE_IN_DIR/acultivar_basnpsmaf.bcf.gz

# Extracting variants and estimating founder and cultivar relatedness.
founder_basnpsmaf_impute_bcf=$PEDIGREE_IN_DIR/founder_basnpsmaf_impute.bcf.gz
founder_pos=$PEDIGREE_IN_DIR/founder_pos.txt
pedig_all_ids=$PEDIGREE_IN_DIR/pedig_all_ids.txt
pedig_ldfilt_bcf=$PEDIGREE_IN_DIR/pedig_ldfilt.bcf.gz
pedig_pos=$PEDIGREE_IN_DIR/pedig_pos.txt
founders_af=$PEDIGREE_IN_DIR/founders_af.txt
founders_winter_af=$PEDIGREE_IN_DIR/founders_winter_af.txt
founders_spring_af=$PEDIGREE_IN_DIR/founders_spring_af.txt
founders_spring_af=$PEDIGREE_IN_DIR/founders_spring_af.txt
founders_all_af=$PEDIGREE_IN_DIR/founders_all_af.txt
pedig_base_af=$PEDIGREE_IN_DIR/pedig_base_af.txt
pedig_ldfilt_plink=$PEDIGREE_IN_DIR/pedig_ldfilt_plink

# Selection scans.
pocultivar_basnpsmaf_impute_vcf=$PEDIGREE_IN_DIR/pocultivar_ldfiltsel.vcf.gz
pocultivar_basnpsmaf_impute_bcf=$PEDIGREE_IN_DIR/pocultivar_basnpsmaf_impute.bcf.gz
pedig_all_basnpsmaf_impute_bcf=$PEDIGREE_IN_DIR/pedig_all_basnpsmaf_impute.bcf.gz

cultivpedig_ba_bcf=$PEDIGREE_IN_DIR/cultivpedig_ba.bcf.gz
founders_ba_bcf=$PEDIGREE_IN_DIR/founders_ba.bcf.gz
founders_winter_af_full=$PEDIGREE_IN_DIR/founders_winter_af_full.txt
founders_semiwinter_af_full=$PEDIGREE_IN_DIR/founders_semiwinter_af_full.txt
founders_spring_af_full=$PEDIGREE_IN_DIR/founders_spring_af_full.txt
founders_all_af_full=$PEDIGREE_IN_DIR/founders_all_af_full.txt
pedig_base_af_full=$PEDIGREE_IN_DIR/pedig_base_af_full.txt
selscan_deltap=$PEDIGREE_IN_DIR/selscan_deltap
selection_windows=$PEDIGREE_OUT_DIR/selection_windows.txt
pedig_ba_bcf=$PEDIGREE_IN_DIR/pedig_ba.bcf.gz
cultivall_ba_bcf=$PEDIGREE_IN_DIR/cultivall_ba.bcf.gz
cultivall_ba_pi="$PEDIGREE_OUT_DIR/cultivall_ba_pi.txt"
cultivall_ba_tajima="$PEDIGREE_IN_DIR/cultivall_ba_tajima"
acultivar_ba_bcf=$PEDIGREE_IN_DIR/acultivar_ba.bcf.gz
acultivar_ba_eupad_bcf=$PEDIGREE_IN_DIR/acultivar_ba_eupad.bcf.gz

# Regions of interest analysis.
a8_basnpsmaf_bcf=$PEDIGREE_IN_DIR/a8_basnpsmaf.bcf.gz
c3_basnpsmaf_bcf=$PEDIGREE_IN_DIR/c3_basnpsmaf.bcf.gz

cactus_fastas="$PEDIGREE_IN_DIR/cactus_fastas/"
merge_filtsub_rename_bcf="$PEDIGREE_IN_DIR/merge_filtsub_rename.bcf.gz"
merge_filtsub_eupad_bcf="$PEDIGREE_IN_DIR/merge_filtsub_eupad.bcf.gz"
eu_founders_bcf=$PEDIGREE_IN_DIR/eu_founders.bcf.gz

rest_ids=$PEDIGREE_IN_DIR/rest_ids.txt
rest_basnpsmaf_bcf=$PEDIGREE_IN_DIR/rest_basnpsmaf.bcf.gz

cultiv_ldfilt_bcf=$PEDIGREE_IN_DIR/cultiv_ldfilt.bcf.gz
pedig_basnpsmaf_bcf=$PEDIGREE_IN_DIR/pedig_basnpsmaf.bcf.gz

pbwt_ref_samples=$PEDIGREE_IN_DIR/pbwtin_ref.txt
pbwt_target_samples=$PEDIGREE_IN_DIR/pbwtin_target.txt
pbwt_all_samples=$PEDIGREE_IN_DIR/pbwtin_all.txt

pbwt_bcf_target_ba=$PEDIGREE_IN_DIR/pbwtin_target_ba.bcf.gz
pbwt_bcf_ref=$PEDIGREE_IN_DIR/pbwtin_ref.bcf.gz
pbwt_bcf_target_sel=$PEDIGREE_IN_DIR/pbwtin_target_sel.bcf.gz
pbwt_bcf_all=$PEDIGREE_IN_DIR/pbwtin_all.bcf.gz

pbwt_positions_target_ba=$PEDIGREE_IN_DIR/pbwtin_positions_target_ba.txt
pbwt_positions_ref=$PEDIGREE_IN_DIR/pbwtin_positions_ref.txt

pbwt_kinship=$PEDIGREE_IN_DIR/pbwt_kinship
