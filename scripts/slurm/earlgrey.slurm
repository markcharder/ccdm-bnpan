#!/bin/bash --login

#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem=30G
#SBATCH --partition=work
#SBATCH --output=logs/eg_%j.out
#SBATCH --error=logs/eg_%j.err
#SBATCH --ntasks=1
#SBATCH --nodes=1

# Use or set scratch environment variable.
MYSCRATCH="${MYSCRATCH:}"
if [ -z "$MYSCRATCH" ]; then
    echo "MYSCRATCH is not set. Set scratch path."
    exit 1
fi

# This may not work with slurm - aim to source scripts/lib/common_functions.sh
SCRIPT_DIR="$(dirname $(dirname "$(realpath "$0")"))"
LIB_DIR="$SCRIPT_DIR"/lib
source "$LIB_DIR"/common_functions.sh
EG_SCRIPT="$SCRIPT_DIR"/earlgrey.slurm

# You can run the script with different options.
sbatch "$EG_SCRIPT" \
    --ram yes \
    --genome-dir $MYSCRATCH/genomes/ \
    --variety mystic \
    --genome-file mystic.fasta \
    --out-dir $MYSCRATCH/earlgrey/
