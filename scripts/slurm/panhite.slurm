#!/bin/bash --login

#SBATCH --job-name=panhite_test
#SBATCH --output=logs/panhite_test.out
#SBATCH --error=logs/panhite_test.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --nodes=1

### Set these variables to run panhite ###

WORKFLOW_DIR="$MYSCRATCH/workflows/"
PANHITE_WORK="$WORKFLOW_DIR/test_panhite/work_tmp/"
PROJECT_DIR="$WORKFLOW_DIR/ccdm-bnpan/"
SOURCE_DIR="$WORKFLOW_DIR/test_panhite/HiTE/"
PAN_GENOMES_DIR="$PROJECT_DIR/data/raw/genomes/"
GENOME_LIST="$PROJECT_DIR/data/conf/genome_list.txt"
OUT_DIR="$PROJECT_DIR/data/processed/panHiTE/"
CONDA_NAME="$WORKFLOW_DIR/test_panhite/HiTE/HiTE_conda"
SCRIPTS_DIR="$PROJECT_DIR/scripts/lib/"

### If this is a test run, point to test data ###

TEST=false
while [[ $# -gt 0 ]]; do

    case $1 in
        --test)
            TEST=true
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
    shift
done

if [ "$TEST" = true ] ; then
    PAN_GENOMES_DIR="$PROJECT_DIR/tests/data/raw/genomes/"
    GENOME_LIST="$PROJECT_DIR/tests/data/conf/genome_list.txt"
    OUT_DIR="$PROJECT_DIR/tests/data/processed/test_panHiTE/"
fi

### Set up nextflow environment ###

# Small function to make sure nextflow is loaded and caches correctly.
source "$SCRIPTS_DIR"/common_functions.sh
nextflow_setup

# Prevent programs from caching in $HOME.
export MPLCONFIGDIR="$MYSCRATCH/.matplotlib_cache"
mkdir -p "$MPLCONFIGDIR"

export CONDA_CACHE="$MYSCRATCH/.conda_cache"
mkdir -p "$CONDA_CACHE"

# Copy custom config to panhite source directory.
cp "$PROJECT_DIR"/conf/panhite.config "$SOURCE_DIR"/nextflow.config

### Run panhite workflow ###

cd $SOURCE_DIR
nextflow \
    run \
        panHiTE.nf \
	    -resume \
            -profile conda,hpc \
            --debug 1 \
            --conda_name "$CONDA_NAME" \
            --conda_cache "$CONDA_CACHE" \
            --pan_genomes_dir "$PAN_GENOMES_DIR" \
            --genome_list "$GENOME_LIST" \
            --out_dir "$OUT_DIR" \
            --work_dir "$PANHITE_WORK"
