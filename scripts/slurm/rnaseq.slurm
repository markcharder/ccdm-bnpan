#!/bin/bash --login

#SBATCH --job-name=rnaseq
#SBATCH --output=logs/rnaseq_%j.out
#SBATCH --error=logs/rnaseq_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --nodes=1

### Set these variables to run panhite ###

WORKFLOW_DIR="$MYSCRATCH/workflows/"
PROJECT_DIR="$WORKFLOW_DIR/ccdm-bnpan/"
SCRIPTS_DIR="$PROJECT_DIR/scripts/lib"

### Set up nextflow environment ###

# Small function to make sure nextflow is loaded and caches correctly.
source "$SCRIPTS_DIR"/common_functions.sh
nextflow_setup
singularity_setup

### Run rnaseq workflow ###

cd "$PROJECT_DIR"
nextflow \
    run \
        . \
        -profile hpc,singularity \
        --workflow "rnaseq"
