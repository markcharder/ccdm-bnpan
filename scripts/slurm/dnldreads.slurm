#!/bin/bash --login

#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/dnldreads_$A_%a.out
#SBATCH --error=logs/dnldreads_%A_%a.err

# Load the aws command line module.
module load awscli/1.29.41

# The file containing the paths of reads should be in $1.
S3_READS_FILE="$1"

# If the file does not exist, exit.
if [[ -z "$S3_READS_FILE" ]]; then
    echo "Provide the file containing paths to reads in the s3 bucket."
    exit 1

fi

# Exit if MYSCRATCH not set.
if [[ -z "$MYSCRATCH" ]]; then
    echo "Set the variable \"$MYSCRATCH\"."
    exit 1

fi

# Set the path for saving the reads files.
WF_PATH="$MYSCRATCH"/workflows/ccdm-bnpan/
READS_PATH="$WF_PATH"/data/graphmap/reads_illumina/
mkdir -p "$READS_PATH"

# Get an array of reads files.
S3_READS_ARR=($(cat "$S3_READS_FILE"))
S3_READS_FULL="${S3_READS_ARR[$SLURM_ARRAY_TASK_ID]}"
S3_READS_BASE="$(basename "$S3_READS_FULL")"

# Run aws s3.
aws s3 cp "$S3_READS_FULL" "$READS_PATH"

# Then touch the file to update its metadata.
# This avoids immediate deletion due to scratch purge policy.
touch "$READS_PATH"/"$S3_READS_BASE"
