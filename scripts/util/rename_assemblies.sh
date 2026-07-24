#!/bin/bash --login

inputs=($(cat "verkko_assem.txt"))
file=${inputs[$SLURM_PROCID]}
vname=$(echo "$file" | cut -d"_" -f 1)
newfile=$(echo ${file/assembly/$vname} | cut -d"/" -f 2)
srun \
    rsync $file $newfile