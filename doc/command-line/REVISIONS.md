# Revisions
These are the steps taken to address reviewers' comments on the manuscript, which
was uploaded as a preprint at: https://doi.org/10.64898/2026.07.29.740455
Unfortunately, the variants used as input were not backed up, so the graph mapping
and variant calling pipeline was re-run. Two additional scripts are included for
downloading the reads from an aws bucket, where they are stored. Note, these reads
are due to be released from NCBI on 1st December under BioProject PRJNA1504441.
These are the scripts run on the HPC to download the reads:
```
scripts/slurm/dnldreads.slurm
scripts/slurm/dnldreads.sub
```
Simply run ```bash scripts/slurm/dnldreads.sub``` to download the reads to
```data/graphmap/reads_illumina```. Then, make the variety IDs config file for the
nextflow pipeline:
```bash
mkdir -p data/graphmap/conf/
python3 scripts/util/make_read_pairs.py
```
Then, run both graphmap and graphcall with:
```bash
sbatch scripts/slurm/graphvar.slurm
```
