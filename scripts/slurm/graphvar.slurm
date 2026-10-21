#!/bin/bash --login

#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --error=logs/graphvar.err
#SBATCH --output=logs/graphvar.out

#####################################################################

# This is the slurm script used to run both graphmap and graphcall.
# nextflow workflows. These take the genomic reads we generated, and
# map them to the pan-genome graph, then call variants from the reads.

#####################################################################

# Load modules.
module load nextflow/25.04.6
module load singularity/4.1.0-slurm

# Export environment variables to control caching.
export NXF_SINGULARITY_CACHEDIR="$MYSCRATCH/.singularity"
export NXF_HOME="$MYSCRATCH/.nextflow"
export SINGULARITY_CACHEDIR="$MYSCRATCH/.singularity"

# Get path of workflow directory.
if [[ ! -s "$MYSCRATCH" ]]; then
    echo "Set the variable \"$MYSCRATCH\"."
    exit 1
fi

# Get to the workflow root path.
export WF_PATH="${1:-"$MYSCRATCH"/workflows/ccdm-bnpan/}"
cd "$WF_PATH"

# Unset sbatch export so child processes inherit the environment.
unset SBATCH_EXPORT

# Check the reads and graph are present.

if [[ ! -s "$WF_PATH"/data/graphmap/graph/ccdm-bnpan.gbz ]]; then
    echo "Graph not present."
    exit 1
fi

if [[ ! -s "$WF_PATH"/data/graphmap/graph/ccdm-bnpan.hapl ]]; then
    echo "Haplotypes not present."
    exit 1
fi

if [[ ! -s "$WF_PATH"/data/graphmap/conf/variety_ids.txt ]]; then
    echo "Variety IDs file not present."
    exit 1
fi

# Also make sure the singularity container is present.
if [[ ! -s "$WF_PATH"/env/containers/graphtools.sif ]]; then
    echo "Container for graphtools not available."
    exit 1
fi

# Run the graphmap pipeline.
nextflow run \
    main.nf \
        -profile singularity,hpc \
        --workflow graphmap

# Run the graphcall pipeline.
nextflow run \
    main.nf \
        -profile singularity,hpc \
        --workflow graphcallg
