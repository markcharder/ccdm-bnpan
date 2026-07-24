# __CCDM _Brassica napus_ pan-genome (ccdm-bnpan)__
This repository contains the code used to annotate a _Brassica napus_ pan-genome assembled at the Centre for Crop and Disease Management (CCDM), Curtin University, Australia. The scripts referenced in the methods section of the manuscript are described in this README.
## Overview
Workflows were developed in three different ways throughout the project, as stand-alone bash scripts, nextflow workflows and markdown notebooks with code snippets.
## Stand-alone bash scripts
### Genome assembly - `ccdm-bnpan/scripts/util/assemble_genomes_execute.sh`
Genomes were assembled on the Setonix supercomputer based in the Pawsey Super-computing Centre in Australia. Scripts for this are in `ccdm-bnpan/scripts/util/`. Three scripts were used to assemble the genomes:
- `assemble_genomes_submit.sh`: submits `assemble_genomes_slurm.sh` to the queue
- `assemble_genomes_slurm.sh`: runs `assemble_genomes_execute.sh`
- `assemble_genomes_execute.sh`: executes the appropriate genome assembly function

Before running `assemble_genomes_submit.sh`, set the variables at the top of the script to point to the reads and output directories and the Verkko conda environment. Reads must be named based on the type of reads. Accepted file names are:
- Nanopore reads = `<VARIETY>_simplex.fastq.gz` and `<VARIETY>_duplex.fastq.gz`
- PacBio reads = `<VARIETY>_hifi.fastq.gz`

Then run:
```bash
bash scripts/util/assemble_genomes_submit.sh
```
### Assembly scaffolding - `ccdm-bnpan/scripts/util/scaffold_genomes.sh`
This scaffolds (without correction) assembled contigs against the reference genome for __B. napus__ ZS11. Set the paths at the top of the file to point to the assemblies from the previous step and the downloaded reference genome. Then run:
```bash
bash scripts/util/scaffold_genomes.sh
```
