# Extract graph regions for A07 sweeps
## Background
On chromosome A07, there are multiple, often quite long, sweeps. A07 seems at first glance to be the most 'selected' chromosome<\br>
Some haplotypes are a bit strange. There are long blocks of SNPs that are all reference whereas the alternate is almost<\br>
all alt alleles for the same block. Weirdly, the all reference blocks appear in some spring lines, like Zephyr, even though<\br>
the reference is an Asian semi-winter and likely more related to other semi-winters(?).<\br>
The objectives here are to:
1) Extract regions surrounding potential sweeps on A07 from the graph GBZ.
2) Use them to assess the allelic variation in these regions.
3) Assess the quality of the variant calls using e.g. read depth to determine whether they are artifacts.
4) Extract gene sequences within them.
5) Using other available genomes, place them in a broader context.
## Code
### Run vg chunk to extract A07 regions

```bash

# Chunk each of the A07 regions only.

# This is to convert accessions to chromosome IDs.
# X5A scaffolds have accession names rather than chromosomes.

# Declare an array linking chromosome to accession.
declare -A CHROM_TO_ACC

# This is a list of chromosome IDs.
chroms=(A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 C01 C02 C03 C04 C05 C06 C07 C08 C09)

# Accessions start at 123876 and increment by 1.
acc=123876
for c in "${chroms[@]}"; do
    CHROM_TO_ACC["$c"]="CP${acc}.1"
    (( acc++ ))
done

# Skip header with tail -n +2
tail -n +2 $selection_windows | while IFS=$'\t' read -r chrom win_start win_end peak_z_af peak_ihs both width_kb source; do

    # Build name to match the graph.
    acc="${CHROM_TO_ACC[$chrom]}"
    path="X5A_ref#0#${acc}"
    slug="${source}_${chrom}_${win_start}_${win_end}"

    # Only do if the chromosome is A07.
    if [ $chrom == A07 ]; then
        echo "Extracting ${path}:${win_start}-${win_end} → ${slug}"

        vg \
            chunk \
                --threads 30 \
                --xg-name "$graph_gbz" \
                --path "${path}:${win_start}-${win_end}" \
                --context-steps 0 \
                --output-fmt gfa \
                > "${HAPLOTYPES_OUT_DIR}/${slug}.gfa"

        START=$(grep "^S" "${HAPLOTYPES_OUT_DIR}/${slug}.gfa" | awk '{print $2}' | sort -n | head -1)
        END=$(grep  "^S" "${HAPLOTYPES_OUT_DIR}/${slug}.gfa" | awk '{print $2}' | sort -n | tail -1)

        vg \
        chunk \
            --threads 30 \
            --xg-name "$graph_gbz" \
            --node-range "${START}:${END}" \
            --output-fmt gfa \
            > "${HAPLOTYPES_OUT_DIR}/${slug}_haplotypes.gfa"

        grep -c "^W\|^P" "${HAPLOTYPES_OUT_DIR}/${slug}_haplotypes.gfa"
    
    fi

done

path='X5A_ref#0#CP123882.1'
win_start=750000
win_end=2000000
slug=test
vg \
    chunk \
        --threads 30 \
        --xg-name "$graph_gbz" \
        --path "${path}:${win_start}-${win_end}" \
        --context-steps 0 \
        --output-fmt gfa \
        > "${HAPLOTYPES_OUT_DIR}/${slug}.gfa"

START=$(grep "^S" "${HAPLOTYPES_OUT_DIR}/${slug}.gfa" | awk '{print $2}' | sort -n | head -1)
END=$(grep  "^S" "${HAPLOTYPES_OUT_DIR}/${slug}.gfa" | awk '{print $2}' | sort -n | tail -1)

vg \
    chunk \
        --threads 30 \
        --xg-name "$graph_gbz" \
        --node-range "${START}:${END}" \
        --output-fmt gfa \
        > "${HAPLOTYPES_OUT_DIR}/${slug}_all.gfa"

```