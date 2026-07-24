
# Set up working environment.
mkdir -p ragtag/{output,log,tmp}
cd ragtag
conda activate ragtag
GENOME_DIR="/path/to/genomes/" # Directory containing assemblies from assemble_genomes_execute.sh
REF_GENOME_DIR="/path/to/reference_genome/" # Directory containing the downloaded reference (ZS11).

# Get a list of final genomes, based on the manual selection process above.
VARIETIES=("Shiralee" "Mutu" "Yickadee" "Drakkar" "Westar" "Maluka" "Marnoo" "Rainbow" "Barossa" "Haya" "Zephyr" "Chisaya" "Eureka" "Castle" "Bronowski" "Dunkeld" "Stubby" "Range" "Insignia" "Trooper" "Chikuzen" "Mystic")
METHODS=("flye" "flye" "verkko" "flye" "flye" "flye" "flye" "flye" "verkko" "verkko" "flye" "verkko" "verkko" "flye" "flye" "verkko" "verkko" "verkko" "verkko" "verkko" "verkko" "verkko")

# Scaffolding using the ZS11 genome.
awk '/^>/ {f=0} /Chromosome/ {f=1} f' \
    "$REF_GENOME_DIR"/GWHANRE00000000.genome.fasta > "$REF_GENOME_DIR"/ZS11_chromosomes.fasta

# Ragtag won't run on gzipped files.
ls "$GENOME_DIR"/*_fasta/*gz | \
    parallel --jobs 32 "gunzip -c {} > tmp/{/.}"

# Then, run ragtag to scaffold the genomes.
for ((i=0;i<${#VARIETIES[@]};i++)); do

    variety=${VARIETIES[$i]}
    method=${METHODS[$i]}

    echo "Running ragtag for $variety and $method on file ../../data/pedigree_genomes/${method}_fasta/${variety}_${method}.fasta.gz"
    ragtag.py \
        scaffold \
        -o output/${variety}_${method}_ragtag \
        -w -u -t 32 \
        "$REF_GENOME_DIR"/ZS11_chromosomes.fasta \
        tmp/${variety}_${method}.fasta > log/${variety}_${method}_ragtag_log.txt 2>&1

done

rm -rf tmp
