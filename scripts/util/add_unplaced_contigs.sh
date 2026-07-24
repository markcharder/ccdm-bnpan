#!/bin/bash

SCRIPTS_DIR="$(dirname "$(dirname "$(realpath "$0")")")"
PROJECT_DIR="$(dirname "$SCRIPTS_DIR")"
FULL_GENOMES="$PROJECT_DIR"/data/raw/full_genomes/
UNPLACED_CONTIGS="$PROJECT_DIR"/data/raw/unplaced_contigs/
mkdir -p "$UNPLACED_CONTIGS"

combine_genomes() {

    for genome in "$FULL_GENOMES"/*; do

        if [[ $(file -b "$genome" | cut -d' ' -f1) == "gzip" ]]; then
            gunzip "$genome"
            local genome="${genome%.*}"
        fi

        genome_base=$(basename "$genome")

        unplaced=$(awk '
            $0 ~ /^>scaffold/ ||
            $0 ~ /^>Scaffold/ ||
            $0 ~ /^>contig/ ||
            $0 ~ /^>tig/ {
                print_it=1
                }
            print_it {
                print
                }
            $0 ~ />/ &&
            !($0 ~ /^>scaffold/ ||
            $0 ~ /^>Scaffold/ ||
            $0 ~ /^>contig/ ||
            $0 ~ /^>tig/) {
                print_it=0
                }' "$genome")

        if [[ -z "$unplaced" ]]; then
            echo "$genome has no unplaced contigs."
        else
            parts=($(echo "${genome_base%.*}" | tr '_' '\n'))

            if [[ ${parts[0]} =~ (bnapus|boleracea|brapa|bnigra) ]]; then
                variety="${parts[1]}"
            else
                variety="${parts[0]}"
            fi

            echo "$unplaced" > "$UNPLACED_CONTIGS"/${variety}_U-genome.fa

        fi

    done

}

combine_genomes