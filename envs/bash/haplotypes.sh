source "$(dirname ${BASH_SOURCE[0]})/pedigree.sh"

HAPLOTYPES_IN_DIR=$WF_DIR/data/haplotypes/inputs/
HAPLOTYPES_OUT_DIR=$WF_DIR/data/haplotypes/out/

mkdir -p $HAPLOTYPES_OUT_DIR
graph_gfa=$HAPLOTYPES_IN_DIR/ccdm-bnpan.gfa
