#!/bin/bash --login

ACCOUNT="${PAWSEY_ACCOUNT:-"none"}"

if [[ "$ACCOUNT" == "none" ]]; then
    echo "No account set. Exiting."
    echo "Set PAWSEY_ACCOUNT env variable."
    exit 1
fi

READS_DIR="/path/to/reads/" # Set to reads directory.
VARIETIES=($(ls $READS_DIR | cut -d"_" -f 1 | sort -u))
N_JOBS=$((${#VARIETIES[@]} - 1))
OUT_DIR="/path/to/output/" # Set to output directory.
TELOMERE="TTTAGGG"
VERKKO_ENV="/path/to/verkko/" # Path to Verkko conda env.

for variety in ${VARIETIES[@]}; do

    assembly_type=""

    if [[ -f $READS_DIR/${variety}_hifi.fastq.gz ]]; then

        if [[ -f $READS_DIR/${variety}_simplex.fastq.gz ]]; then
            assembly_type="pbnp"
        else
            assembly_type="pb"
        fi

    else
        assembly_type="np"

    fi

    echo "$variety $assembly_type"

done > assembly_types.txt

# Submit the job to slurm.
mkdir -p log
sbatch \
    --array=0-$N_JOBS \
    --partition=work \
    --account="$ACCOUNT" \
    --exclusive \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=64 \
    --job-name=assemble_genomes \
    --output=log/assemble_genomes_%A_%a.out \
    --error=log/assemble_genomes_%A_%a.err \
    --time=24:00:00 \
    assemble_genomes_slurm.sh \
    assembly_types.txt \
    "$READS_DIR" \
    "$OUT_DIR" \
    "$TELOMERE" \
    "$VERKKO_ENV"
