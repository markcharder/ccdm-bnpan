# Dependencies
2.28-r1209

# Environment
```bash
source envs/bash/genetests.sh
```

# Code
## Transfer annotations from ragtag scaffolds to contigs
When braker3 was run the first time around, it ran against genomes scaffolded to zs11. Since the graph contains<br>
unscaffolded assemblies in their constituent contigs, I need to transfer the annotations from braker3 (and panhite)<br>
to the contigs to find their context in the pan-genome graph.<br>
I have decided to align contigs to the zs11 scaffolds for each genome and find their exact matches. Based on the<br>
coordinates of exact matches, I should be able to change the coordinates of the annotations.
### Unpack annotations (if not done already)
The annotations were generated on the HPC cluster. I will unpack them in the data directory.<br>
Move on if this was done already.<br>
```bash

cd $DATA/annotations/genes/
tar xvf $DATA/workflows/outputs/backup/2025-12-22_braker3.tar.gz
braker_path=`pwd`/braker3
cd $DATA/workflows/outputs/backup
tar xvf $DATA/workflows/outputs/backup/2026-01-27_cactus-in.tar.gz
cactus_path=`pwd`/data

wfcd
cd ccdm-bnpan/data/genetests/inputs

# This is to get the outputs of braker3.
for dir in `ls -d $braker_path/*/*`; do

    if [[ $dir != *"augustus"* ]]; then
        echo "Processing $dir."
        suffix=${dir##*/}
        rsync --progress $dir/braker.gtf ${suffix}_braker3.gtf
        rsync --progress $dir/braker.codingseq ${suffix}_braker3-cds.fasta
        rsync --progress $dir/braker.aa ${suffix}_braker3-aa.fasta
    fi

done

# And the inputs of cactus.
ln -s $cactus_path/*gz .

# These should match the inputs to braker3 / panhite but they are not soft-masked.
ln -s $DATA/assemblies/scaffold_renamed/non-panSN/* .

```
### Align contigs to scaffolds
Then use minimap2 to align cactus contigs to scaffolded genomes.<br>
```bash
CONTIG_NAMES=($(ls data/genetests/inputs/*gz))

minimap_fun(){

    contigs=$1
    data_dir=$2
    out_dir=$3
    threads=$4
    ram=$5

    file_base=$(basename $contigs)
    variety_name=$(echo $file_base | cut -f 1 -d "_")
    scaffolds=($(ls $data_dir/${variety_name}_genome-complete.fa 2> /dev/null))

    if [[ -f "${scaffolds[0]}" ]]; then
        echo -e "Scaffolds: $scaffolds Variety: $variety_name Contigs: $contigs" >&2
        minimap2 \
            -o $out_dir/${variety_name}_mm2.paf \
            -t $threads \
            -x asm5 \
            --cs \
            -c \
            $scaffolds $contigs

    else
        echo -e "Variety: $variety_name has no scaffolds" >&2
    fi
}

export -f minimap_fun

mkdir -p $GENETESTS_OUT_DIR/minimap2
mkdir -p logs

printf "%s\n" ${CONTIG_NAMES[@]} | \
    parallel --jobs 2 "minimap_fun {} data/genetests/inputs data/genetests/out/minimap2 15 1> logs/{/.}_mm2_log.txt 2>&1"

```

### Transfer the annotations
```bash

mkdir -p $GENETESTS_OUT_DIR/annotations
PAFS=($(ls $GENETESTS_OUT_DIR/minimap2/*paf))

for paf in ${PAFS[@]}; do
    prefix=$(basename $paf | cut -d"_" -f 1)
    gtf=$GENETESTS_IN_DIR/${prefix}_braker3_braker3.gtf
    echo $paf $prefix $gtf
    contigs=$(ls $GENETESTS_IN_DIR/${prefix}*gz)

    if [ -f $contigs ]; then
        
        if [ ! -f ${contigs%.*} ]; then
            zcat $contigs > ${contigs%.*}

        fi

        ./util/transfer_gtf_scaffold_to_contigs.py \
            --paf $paf \
            --gtf $gtf \
            --output $GENETESTS_OUT_DIR/annotations/$prefix \
            --test \
            --contigs ${contigs%.*} \
            --cds $GENETESTS_IN_DIR/${prefix}_braker3_braker3-cds.fasta \
            1> logs/${prefix}_annot-trans.txt 2>&1

    else
        echo "Contigs not found for $prefix"
    fi

done

```

### Extract genes from selected regions
#### With vg rna
Option 1 is to add transcript paths to the graph to find all genes overlapping swept windows<br>
This is not straightforward. It could be revisited at a later date.<br>
```bash

# Find start and end nodes for selected regions.
declare -A CHROMOSOME_IDS

while IFS=$'\t' read -r acc rid fid bcfid length; do
    CHROMOSOME_IDS["$rid"]="X5A_ref#0#$acc"
done < doc/info/chromosomes_rename_key_lens.txt

mkdir -p $GENETESTS_OUT_DIR/selection_windows/
while IFS=$'\t' read -r chrom winstart winend rest; do

    if [ "$chrom" != "CHROM" ]; then
        acc=${CHROMOSOME_IDS["$chrom"]}
        echo -e "$acc\t$((winstart - 1))\t$winend"            

    fi

done < data/pedigree/out/selection_candidates_merged.txt | sort -V -k1,1 -k2,2 > $GENETESTS_OUT_DIR/selection_windows/selection_candidates_merged.bed

# Save to subgraphs.
vg \
    find \
        --path-bed $GENETESTS_OUT_DIR/selection_windows/selection_candidates_merged.bed \
        --save-to $GENETESTS_OUT_DIR/selection_windows/selwindows_ \
        --xg-name $graph_gbz

```
#### With minimap2 alignments
An alternative option is to align all contigs for one arbitrarily good reference genome<br>
to X5A chromosomes. Then, any transcripts within unique long-range alignments on X5A that<br>
are within selected windows can be considered.<br>
```bash
# Extract chromosome sub-problems.
CHROMSP=($(tar -tf $DATA/workflows/outputs/backup/2026-01-23_14-04-00_cactus.tar.gz | grep -E "Chikuzen|X5A" | grep -v AMBIGUOUS))
tar xvf $DATA/workflows/outputs/backup/2026-01-23_14-04-00_cactus.tar.gz "${CHROMSP[@]}"
ACCS=($(printf "%s\n" ${CHROMSP[@]} | cut -d"/" -f 3 | sort -u))
for acc in ${ACCS[@]}; do

    minimap2 \
        -t 12 \
        -x asm5 \
        --cs \
        "2026-01-23_14-04-00_cactus/chrom-subproblems/$acc/fasta/X5A_ref_$acc.fa.gz" \
        "2026-01-23_14-04-00_cactus/chrom-subproblems/$acc/fasta/Chikuzen_verkko_$acc.fa.gz" \
        > $GENETESTS_OUT_DIR/minimap2/chik2x5a_$acc.paf

done
```