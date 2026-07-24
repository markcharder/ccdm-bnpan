#!/bin/bash --login

### Get input parameters ###
mapfile -t ASSEMBLY_COMMANDS < $1
read -r -a CURRENT_COMMAND <<< "${ASSEMBLY_COMMANDS[$SLURM_ARRAY_TASK_ID]}"

READS_DIR=$2
OUT_DIR=$3
TELOMERE=$4
VERKKO_ENV=$5

### Set up working environment ###
conda activate "$VERKKO_ENV"

### Pass parameters to execute script ###
srun \
    --export=ALL \
    bash \
        assemble_genomes_execute.sh \
        "${CURRENT_COMMAND[0]}" \
        "${CURRENT_COMMAND[1]}" \
        "$READS_DIR" \
        "$OUT_DIR" \
        "$TELOMERE"
