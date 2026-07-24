#!/usr/bin/env bash

OUT_FILE="chromosome_ids.txt"
DELIM="_"

if [ $# -eq 0 ]; then
    echo "Usage: chrominfo.sh --fasta_dir <directory> [--out_file <file> --delim <delimiter>]"
    exit 1
fi

while [ $# -gt 0 ]; do

    case "$1" in

        --fasta_dir)
            FASTA_DIR="$2"
            shift 2
            ;;
        
        --delim)
            DELIM="$2"
            shift 2
            ;;

        --out_file)
            OUT_FILE="$2"
            shift 2
            ;;

        --help)
            echo "Usage: chrominfo.sh --fasta_dir <directory> [--out_file <file> --delim <delimiter>]"
            exit 0
            ;;

        *)
            echo "Unknown argument: $1"
            exit 1
            ;;

    esac

done

if [ ! -d "$FASTA_DIR" ]; then
    echo "Error: Directory $FASTA_DIR does not exist." >&2
    exit 1
fi

if [ -f "$OUT_FILE" ]; then
    echo "Warning: outfile already exists, overwriting." >&2
    rm "$OUT_FILE"
fi

FASTA_FILES=($(ls $FASTA_DIR/*))

if [ ${#FASTA_FILES[@]} -eq 0 ]; then
    echo "Error: No fasta files found in directory \"$FASTA_DIR\"." >&2
    exit 1
fi

CHROM_IDS=($(ls $FASTA_DIR/ | cut -d"$DELIM" -f 1))

for ((i=0;i<${#FASTA_FILES[@]};i++)); do

    file=$(basename "${FASTA_FILES[$i]}")
    chrom_id="${CHROM_IDS[$i]}"

    if [[ "$file" == *.gz ]]; then
        nhaps="$(zcat "$FASTA_DIR/$file" | grep -c "^>")"
    
    else
        nhaps="$(grep -c "^>" "$FASTA_DIR/$file")"
    
    fi

    echo -e "$chrom_id\t$file\t$nhaps" >> "$OUT_FILE"
done