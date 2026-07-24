#!/bin/bash --login

#SBATCH --time=05:00:00
#SBATCH --array=0-397
#SBATCH --output=vcfbub.output
#SBATCH --error=vcfbub.err
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1

module load singularity/4.1.0-nompi

OUT_PATH=/path/to/output # Path for output files

mkdir -p ${OUT_PATH}

bash vcfbub.run $OUT_PATH
