#!/bin/bash --login

#SBATCH --time=05:00:00
#SBATCH --array=0-397
#SBATCH --output=bgzip.output
#SBATCH --error=bgzip.err
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load singularity/4.1.0-nompi
unset SBATCH_EXPORT
export SINGULARITY_CACHEDIR="$MYSCRATCH/.singularity"
# Run this before to improve efficiency:
# singularity pull "docker://stephenturner/bgzip:latest"

bash bgzip.run
