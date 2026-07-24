#!/bin/bash --login

### Get input parameters ###
VARIETY=$1
ASSEMBLY_TYPE=$2
READS_DIR=$3
OUT_DIR=$4
TELOMERE=$5

PBREADS=$READS_DIR/${VARIETY}_hifi.fastq.gz
NP_HERRO=$READS_DIR/${VARIETY}_corrected.fasta.gz
NP_DUPLEX=$READS_DIR/${VARIETY}_duplex.fastq.gz
NP_SIMPLEX=$READS_DIR/${VARIETY}_simplex.fastq.gz

echo "Running pipeline for $VARIETY with assembly type ${ASSEMBLY_TYPE}."

### Set up environment ###
TMPDIR="$(pwd)/tmp_${SLURM_JOB_ID}_$SLURM_ARRAY_TASK_ID"
mkdir -p "$TMPDIR"

SNAKEMAKE_OUTPUT_CACHE="$(pwd)/tmp_sm_${SLURM_JOB_ID}_$SLURM_ARRAY_TASK_ID"
mkdir -p "$SNAKEMAKE_OUTPUT_CACHE"

XDG_CACHE_HOME="$(pwd)/tmp_xdg_${SLURM_JOB_ID}_$SLURM_ARRAY_TASK_ID"
mkdir -p "$XDG_CACHE_HOME"

export TMPDIR
export SNAKEMAKE_OUTPUT_CACHE
export XDG_CACHE_HOME

unset PYTHONHOME
unset PYTHONPATH

### Define assembly functions ###
verkko_pbnp(){

    # Include nanopore and pacbio.
    verkko \
        -d "$READS_DIR/${VARIETY}_nppb" \
        --hifi "$PBREADS" "$NP_HERRO" "$NP_DUPLEX" \
        --nano "$NP_SIMPLEX" \
        --telomere-motif "$TELOMERE" \
        --local-cpus "$SLURM_CPUS_PER_TASK" \
        --haploid
  
}

verkko_pb(){

    # Include pacbio only.
    verkko \
        -d "$READS_DIR/${VARIETY}_pb" \
        --hifi "$PBREADS" \
        --no-nano \
        --telomere-motif "$TELOMERE" \
        --local-cpus "$SLURM_CPUS_PER_TASK" \
        --haploid

}

verkko_np() {

    # Include nanopore only.
    verkko \
        -d "$READS_DIR/${VARIETY}_np" \
        --hifi "$NP_HERRO" "$NP_DUPLEX" \
        --nano "$NP_SIMPLEX" \
        --telomere-motif "$TELOMERE" \
        --local-cpus "$SLURM_CPUS_PER_TASK" \
        --haploid

}

# Check assembly type and run accordingly.
if [ "$ASSEMBLY_TYPE" == "pbnp" ]; then
    verkko_pbnp

elif [ "$ASSEMBLY_TYPE" == "pb" ]; then
    verkko_pb

elif [ "$ASSEMBLY_TYPE" == "np" ]; then
    verkko_np

else
    echo "Unknown assembly type: $ASSEMBLY_TYPE"
    exit 1

fi

rm -rf "$TMPDIR" "$SNAKEMAKE_OUTPUT_CACHE" "$XDG_CACHE_HOME"