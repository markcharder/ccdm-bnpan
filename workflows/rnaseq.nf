#!/usr/bin/env nextflow

def rnaseqHelpMessage() {

    log.info"""

    RNASEQ WORKFLOW
    ============================
    This workflow performs RNA-Seq read alignment to genomes using HISAT2 or downloads and aligns optimal RNA-Seq data
    for using VARUS.

    Optional parameters:


        GENERAL:

        --help                  Show this help message and exit.

        --aligner               Aligner to use: 'hisat2' or 'varus'. Default: 'varus'.
                                If 'varus' is selected, RNA-Seq data will be downloaded automatically from SRA.
                                If 'hisat2' is selected, RNA-Seq reads in --rna_seq_dir are used.

        --rna_seq_dir           Directory containing RNA-Seq reads in FASTQ format (only used for hisat2).
                                Default: ./data/raw/rnaseq/

        --genomes_dir           Directory containing the genomes to annotate.
                                Default: ./data/raw/genomes/

        --repeats_gffs          GFF files containing repeats for masking genomes.
                                Default: "./data/processed/panHiTE/pan_annotate_genomes/*/*{fa,fasta}.sorted.gff"

        --use_unplaced          Whether or not to incorporate unplaced contigs into the genomes.
                                Default: true

        --unplaced_contigs_dir  Directory containing unplaced contigs for merging with placed contigs.
                                Default: "./data/raw/unplaced_contigs/"

        --mask_unplaced         Whether or not to mask unplaced contigs before incorporating.
                                Default: true

        --te_library            Directory containing the transposable element library to mask unplaced contigs.
                                Default: "./data/processed/panTE/panTE.fa"

        --rnaseq_out_dir        Output directory for RNA-Seq workflow results.
                                Default: ./results/rnaseq/

        --index_only            Whether to stop after building the hisat2 index.
                                Default: false

        VARUS:

        --varus_genus           Genus name for VARUS.
                                Default: Brassica

        --varus_vsk             Variety -> species key for VARUS. Optional comma-separated
                                list of mappings between varieties and species names.
                                Example: 'PM:juncea,Westar:napus'. Overrides --varus_species if
                                set.
                                Default: 'PM:juncea,CN115125:nigra,W03:oleracea,Chiifu:rapa,ECD04:rapa'

        --varus_species         Species name for VARUS.
                                Default: 'napus'
        
        --varus_batchSize       Batch size for VARUS.
                                Default: 50000

        --varus_maxBatches      Maximum number of batches for VARUS.
                                Default: 750

        --varus_use_tmp         Whether or not to run VARUS in /tmp. On some systems, this may mean
                                using RAM in a tmpfs backed file system, which could speed up IO.
                                Default: true

""".stripIndent()

}

process filter_small_contigs {

    label 'small'
    tag "${variety}"

    input:
    tuple val(variety), path(unplaced_fasta)

    output:
    tuple val(variety), path("${unplaced_fasta.baseName}-filtered.fa"), emit: filtered_unplaced_fasta

    script:
    """
    remove_small_contigs.py \
        --fasta ${unplaced_fasta} \
        --min-size 1000 > ${unplaced_fasta.baseName}-filtered.fa
    """

    stub:
    """
    touch "${unplaced_fasta.baseName}-filtered.fa"
    """

}

process mask_unplaced_contigs {

    label 'repeatmasker'
    tag "${variety}"

    storeDir "${params.rnaseq_out_dir}/masked_unplaced_contigs/"

    input:
    tuple val(variety), path(genome_fasta), path(repeat_library)

    output:
    tuple val(variety), path("${genome_fasta}.masked"), path("${genome_fasta}.out.gff"), emit: masked_unplaced_contigs

    script:
    """
    forks=\$((${task.cpus} / 4))
    RepeatMasker \
        -parallel \$forks \
        -s \
        -lib ${repeat_library} \
        -nocut \
        -noisy \
        -dir . \
        -xsmall \
        -gff \
        ${genome_fasta}
    """

    stub:
    """
    touch "${genome_fasta}.masked" "${genome_fasta}.out.gff"
    """

}

process combine_variety_files {

    label 'small'
    tag "${variety}"
    
    input:
    tuple val(variety), path(gff_files), path(fasta_files)
    
    output:
    tuple val(variety), path("${variety}_genomes-combined.fa"), path("${variety}_repeats-combined.gff"), emit: cat_files
    
    script:
    """
    # Concatenate all FASTA files for this variety
    cat ${fasta_files} > ${variety}_genomes-combined.fa
    
    # Concatenate all GFF files for this variety
    cat ${gff_files} > ${variety}_repeats-combined.gff
    """

    stub:
    """
    touch "${variety}_genomes-combined.fa"
    touch "${variety}_repeats-combined.gff"
    """
}

process combine_placed_and_unplaced {

    label 'small'
    tag "${variety}"

    storeDir "${params.rnaseq_out_dir}/masked_genomes_full/"

    input:
    tuple val(variety), path(fasta_files)

    output:
    tuple val(variety), path("${variety}_full-genome.fa"), emit: cat_files

    script:
    log.info("Got input: ${variety}, ${fasta_files}")
    """
    # Concatenate all FASTA files for this variety
    cat ${fasta_files} > ${variety}_full-genome.fa
    """

    stub:
    """
    touch "${variety}_full-genome.fa"
    """

}

process softmask_genome {

    label 'small'
    label 'bedtools'
    tag "${variety}"

    storeDir "${params.rnaseq_out_dir}/masked_genomes"

    input:
    tuple val(variety), path(genome_fasta), path(repeats_gff)

    output:
    tuple val(variety), path("${variety}.masked.fa"), emit: masked_genome

    script:
    log.info("Process softmask_genome: storeDir = ${params.rnaseq_out_dir}/masked_genomes/")
    """
    bedtools maskfasta -soft -fi "${genome_fasta}" -bed "${repeats_gff}" -fo "${variety}.masked.fa"
    """

    stub:
    """
    touch "${variety}.masked.fa"
    """

}

process hisat2_index{

    label 'hisat2'
    tag "${variety}"

    // Somewhat resource-intensive step, so keep it.
    storeDir "${params.rnaseq_out_dir}/hisat2_index/"

    input:
    tuple val(variety), path(masked_genome)

    output:
    tuple val(variety), path("${variety}_h2idx.*"), emit: h2idx

    script:
 
    log.info("Building hisat2 index for variety ${variety} with genome ${masked_genome}.")
    """
    hisat2-build \
        --threads ${task.cpus} \
        "${masked_genome}" \
        "${variety}_h2idx"
    """

    stub:
    """
    touch "${variety}_h2idx.1.ht2"
    touch "${variety}_h2idx.2.ht2"
    touch "${variety}_h2idx.3.ht2"
    """
}

process hisat2_align {

    label 'hisat2'
    tag "${variety}-${sample_id}"

    input:
    tuple val(variety), path(h2idx), val(sample_id), path(rna_seq_reads)

    output:
    tuple val(variety),
          val(sample_id),
          path("${variety}-${sample_id}.sam"),
          emit: alignment
    
    script:
    def is_paired = rna_seq_reads instanceof List && rna_seq_reads.size() == 2
    def h2idx_base = "${variety}_h2idx"
    
    if (is_paired) {
        log.info("Aligning paired-end reads for sample ${sample_id} on variety ${variety}.")
        """
        hisat2 \
            --threads ${task.cpus} \
            --dta \
            -x ${h2idx_base} \
            -1 ${rna_seq_reads[0]} \
            -2 ${rna_seq_reads[1]} \
            -S ${variety}-${sample_id}.sam
    """
    } else {
        log.info("Aligning single-end reads for sample ${sample_id} on variety ${variety}.")
        """
        hisat2 \
            --threads ${task.cpus} \
            --dta \
            -x ${h2idx_base} \
            -U ${rna_seq_reads} \
            --sam ${variety}-${sample_id}.sam
    """
    }

    stub:
    """
    touch "${variety}-${sample_id}.sam"
    """
}

process samtools_sort {

    label 'samtools'
    tag "${variety}-${sample_id}"

    // Store it because we may need again.
    storeDir "${params.rnaseq_out_dir}/alignments/"

    input:
    tuple val(variety), val(sample_id), path(sam_file)

    output:
    tuple val(variety),
          val(sample_id),
          path("${variety}-${sample_id}_h2.sorted.bam"),
          emit: sorted_bam

    script:
    """
    samtools view -b -@ ${task.cpus} ${sam_file} | \
    samtools sort -o ${variety}-${sample_id}_h2.sorted.bam
    """

    stub:
    """
    touch "${variety}-${sample_id}_h2.sorted.bam"
    """

}

process varus {

    label 'varus'
    tag "${variety}"

    storeDir "${params.rnaseq_out_dir}/alignments/"

    input:
    tuple val(variety),
        path(masked_genome),
        path(h2idx),
        val(variety_species_map)

    output:
    tuple val(variety), path("${variety}_varus.bam"), emit: varus_bam

    script:
    log.info("Running varus with ${task.cpus} threads and ${task.memory} memory.")
    def species = variety_species_map && variety_species_map.containsKey(variety) ?
                  variety_species_map[variety] : 
                  "${params.varus_species}"
    def genus = "${params.varus_genus}"
    def sp = species

    log.info("Found species for variety ${variety}: ${genus} ${sp}.")

    """

    # Get process nextflow dir.
    NFDIR="\$(pwd)"

    if [ "${params.varus_use_tmp}" = "true" ]; then

        # Set working directory and create it.
        WORKDIR="/tmp/varus_${variety}/"
        mkdir -p \$WORKDIR

        # Copy the masked genome to /tmp.
        cp "${masked_genome}" \$WORKDIR/

    else

        WORKDIR="\$NFDIR"

    fi

    # Make a directory for the varus index.
    mkdir -p \$WORKDIR/${genus}_${sp}/genome/

    # Copy hisat2 index files to the directory from previous run.
    echo "Copying ${h2idx}"
    for file in ${h2idx}; do
        ext="\${file##*_h2idx.}"
        cp "\$file" "\$WORKDIR/${genus}_${sp}/genome/hisatidx.\$ext"
    done

    # Change to working directory.
    cd \$WORKDIR

    # Get species parameters for varus.
    varus_genus=${genus}
    varus_sp=${sp}

    # Set up the varus parameters file.
    echo "--batchSize ${params.varus_batchSize}" > VARUSparameters.txt
    echo "--blockSize 5000" >> VARUSparameters.txt
    echo "--components 1" >> VARUSparameters.txt
    echo "--cost 0.001" >> VARUSparameters.txt
    echo "--deleteLater 0" >> VARUSparameters.txt
    echo "--estimator 2" >> VARUSparameters.txt
    echo "--exportObservationsToFile 1" >> VARUSparameters.txt
    echo "--exportParametersToFile 1" >> VARUSparameters.txt
    echo "--fastqDumpCall fastq-dump" >> VARUSparameters.txt
    echo "--genomeDir ./genome/" >> VARUSparameters.txt
    echo "--lambda 10.0" >> VARUSparameters.txt
    echo "--lessInfo 1" >> VARUSparameters.txt
    echo "--loadAllOnce 0" >> VARUSparameters.txt
    echo "--maxBatches ${params.varus_maxBatches}" >> VARUSparameters.txt
    echo "--mergeThreshold 10" >> VARUSparameters.txt
    echo "--outFileNamePrefix ./${variety}_varus" >> VARUSparameters.txt
    echo "--pathToParameters ./VARUSparameters.txt" >> VARUSparameters.txt
    echo "--pathToRuns ./" >> VARUSparameters.txt
    echo "--pathToVARUS /opt/VARUS/Implementation" >> VARUSparameters.txt
    echo "--profitCondition 0" >> VARUSparameters.txt
    echo "--pseudoCount 1" >> VARUSparameters.txt
    echo "--qualityThreshold 5" >> VARUSparameters.txt
    echo "--randomSeed 1" >> VARUSparameters.txt
    echo "--readParametersFromFile 1" >> VARUSparameters.txt
    echo "--runThreadN ${task.cpus}" >> VARUSparameters.txt
    echo "--verbosityDebug 1" >> VARUSparameters.txt

    # Run VARUS.
    /opt/VARUS/runVARUS.pl \
        --aligner=HISAT \
        --readFromTable=0 \
        --createindex=0 \
        --latinGenus=\$varus_genus \
        --latinSpecies=\$varus_sp \
        --speciesGenome=${masked_genome.name} \
        --logfile=${variety}_varus.log 1> \
        \$NFDIR/${variety}_varus.out 2> \
        \$NFDIR/${variety}_varus.err

    mv \$WORKDIR/${genus}_${sp}/VARUS.bam \$NFDIR/${variety}_varus.bam

    if [ "${params.varus_use_tmp}" = "true" ]; then

        # Clean up /tmp and keep log file.
        mv ${variety}_varus.log \$NFDIR/
        rm -rf \$WORKDIR

    fi
    """

    stub:
    """
    touch "${variety}_varus.bam"
    """

}

workflow RNASEQ {

    // If the help paramter is set, print help for this workflow then exit.
    if (params.help) {
        rnaseqHelpMessage()
        System.exit(0)
    }

    // Take the GFF files from panHiTE annotations and extract the variety name,
    // which is the first part of the '_' separated file name.
    channel
        .fromPath(
            "${params.repeats_gffs}",
                checkIfExists: true
                )
            .map {
                gff ->
                def variety = gff.baseName.split("_")[0]
                tuple(variety, gff)
            }
            .groupTuple()
            .set { repeats_by_variety } // Each variety maps to a list of gff files for each subgenome.

    channel
        .fromPath(
            "${params.genomes_dir}/*.fa",
            checkIfExists: true
            )
        .map { file ->
            def variety = file.baseName.split('_')[0]
            tuple(variety, file)
        }
        .groupTuple()
        .set { genomes_by_variety } // Each variety maps to a list of fasta files for each subgenome.

    // This should make a tuple with each variety pointing to a list of gffs and fastas for subgenomes.
    repeats_by_variety
        .join(
            genomes_by_variety
            )
        .set { combined_varieties }

    // This concatenates gffs and fastas into a single file for
    // each variety, with the suffixes '_genomes-combined.fa' and
    // '_repeats-combined.gff'.
    combine_variety_files(combined_varieties)

    // Create a new directory for caching masked genomes.
    log.info("Using storeDir for masked genomes: ${params.rnaseq_out_dir}/masked_genomes/")
    new File(
        "${params.rnaseq_out_dir}/masked_genomes"
        )
        .mkdirs()
    
    // Softmask the genomes.
    // Takes variety -> gff, fasta
    // emits variety -> masked_fasta
    softmask_genome(
        combine_variety_files.out.cat_files
        )

    // If we want to use unplaced contigs.
    if (params.use_unplaced) {

        // Find them from the input path.
        channel
            .fromPath(
                "${params.unplaced_contigs_dir}/*.{fa,fasta}",
                checkIfExists: true
                )
            .map {
                file ->
                def variety = file.baseName.split('_')[0]
                tuple(variety, file)
            }
            .set { filter_unplaced_input_ch } // This should contain variety -> unplaced contigs.

        // Make sure they have been found.
        log.info("Found unplaced contig files: ${filter_unplaced_input_ch.collect()}")

        // This should remove all contigs below a certain threshold (1000 bp).
        filter_small_contigs(filter_unplaced_input_ch) // emits variety -> filtered_fasta (suffix -filtered.fa)

        // So we have filtered the unplaced contigs.
        // Now, if they need to be masked, run this.
        if (params.mask_unplaced) {

            // Get the repeat library (panTE.fa).
            channel
                .fromPath("${params.te_library}", checkIfExists: true)
                .set { repeat_library_ch }

            log.info("Masking filtered contigs with repeat library: ${repeat_library_ch.collect()}")

            // Take the variety -> filtered_fasta emission from filter_small_contigs.
            filter_small_contigs.out.filtered_unplaced_fasta
                .combine(repeat_library_ch) // Every variety should now map to filtered_fasta and repeat_library.
                .set { mask_unplaced_input_ch }

            // Mask the unplaced contigs with repeatmasker.
            // This emits variety -> masked_fasta, gff (genome_file.fa.masked and genome_file.fa.out.gff)
            mask_unplaced_contigs(mask_unplaced_input_ch)

            // Map the output to remove the gff.
            mask_unplaced_mapped = mask_unplaced_contigs.out.masked_unplaced_contigs
                .map {
                    variety, masked_fasta, _gff ->
                    tuple(variety, masked_fasta)
                }

            // Take the softmasked genome output and perform a left join with mask_unplaced_mapped.
            softmask_genome.out.masked_genome
                .join(
                    mask_unplaced_mapped,
                    remainder: true
                )
                .map{
                    variety, placed, unplaced ->
                    // If there are no unplaced contigs, just return the placed contigs.
                    if (unplaced == null) {
                        tuple(variety, [placed])

                    } else {
                        tuple(variety, [placed, unplaced])
                    }
                }
                .set { full_genomes_ch }

            // Then concatenate the placed and unplaced masked contigs.                
            combine_placed_and_unplaced(full_genomes_ch)
                .set { index_input_ch }

            log.info("Masked and joined unplaced contigs with masked placed contigs.")

        } else {

            // If not masking, just take the filtered unplaced contigs.
            softmask_genome.out.masked_genome
                .join(
                    filter_small_contigs.out.filtered_unplaced_fasta,
                    remainder: true
                )
                .map {
                    variety, placed, unplaced ->
                    // If unplaced is null, just return placed genome
                    if (unplaced == null) {
                        tuple(variety, [placed])
                    } else {
                    // If unplaced exists, combine both
                        tuple(variety, [placed, unplaced])
                    }
                }
                .set { full_genomes_ch }

            // Then combine the masked placed and unmasked unplaced contigs.
            combine_placed_and_unplaced(full_genomes_ch)
                .set { index_input_ch }
            
            log.info("Joined unasked unplaced contigs with masked placed contigs.")

        }

    // If not using unplaced contigs, just use the softmasked genomes.
    } else {

        log.info("Not using unplaced contigs for RNA-Seq alignment.")
        softmask_genome.out.masked_genome
            .map {
                variety, masked_fasta ->
                tuple(variety, masked_fasta)
            }
            .set { index_input_ch }
    }

    log.info("Using storeDir for hisat2 index: ${params.rnaseq_out_dir}/hisat2_index/")

    // Create the cache directory for the hisat2 index.
    new File(
        "${params.rnaseq_out_dir}/hisat2_index"
        )
        .mkdirs()

    // Index the genomes with hisat2.
    // Emits h2idx, which is variety -> hisat2 index files.
    hisat2_index(
        index_input_ch
        )

    // Make the cache directory for alignments in BAM format.
    log.info("Using storeDir for alignments: ${params.rnaseq_out_dir}/alignments/")
    new File(
        "${params.rnaseq_out_dir}/alignments"
        )
        .mkdirs()

    // If we are aligning as well as making just the index.
    if (!params.index_only) {

        // If the aligner is just hisat2.
        if (params.aligner == 'hisat2') {

            // The rnaseq input files need to exist.
            if (!file(params.rna_seq_dir).exists()) {
                error "RNA-Seq directory specified with --rna_seq_dir does not exist: ${params.rna_seq_dir}"
            }

            // Assuming they exist, find paired and single-end reads.
            channel
                .fromFilePairs(
                    "${params.rna_seq_dir}/*_{,R}{1,2}.{fastq,fq}{,gz}"
                    )
                .set { paired_reads_ch }    

            channel
                .fromPath(
                    "${params.rna_seq_dir}/*.{fastq,fq}{,gz}"
                    )
                .filter {
                    file -> 
                    !file.name.matches(/.*_R?[12]\..*/)
                }
                .map {
                    file ->
                    def sample_id = file.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
                    [sample_id, file]
                }
                .set { single_reads_ch }

            // Mix paired and single reads.
            paired_reads_ch
                .mix(single_reads_ch)
                .set { all_reads_ch }

            // Combine with hisat2 index.
            hisat2_index.out.h2idx
                .combine(
                    all_reads_ch
                    )
                .set { h2idx_reads_ch }

            // Align with hisat2.
            hisat2_align(h2idx_reads_ch)

            // Convert to bam and sort with samtools.
            samtools_sort(hisat2_align.out.alignment)

        } else if (params.aligner == 'varus') {

            // If using varus rather than hisat2 directly.
            // Join genome_fasta files with their indexes to make variety -> genome_fasta, index files.
            index_input_ch
                .join(
                    hisat2_index.out.h2idx
                )
                .set { genome_index_ch }

            // Get a mapping between variety and varus species.
            // The only non-napus species is PM, which is a B. juncea genome.
            def variety_species_map = [:]

            // If this option is used (by default it always will be, actually).
            if (params.varus_vsk) {

                // Split the comma-separated list into : separated pairs.
                params.varus_vsk.split(',')
                .each {
                    pair ->
                    def (variety, species) = pair.split(':')
                    variety_species_map[variety.trim()] = species.trim()
                } // Split the pairs to make a dictionary of variety -> species.

            }

            // Take the variety -> genome_fasta, hisat2_index files channel.
            genome_index_ch
                .map {
                    variety, masked_genome, h2idx ->
                    tuple(variety, masked_genome, h2idx, variety_species_map)
                } // Map this to the species map.
                .set { varus_input_ch }
                varus(varus_input_ch) // Run varus on this channel.

        }

    }

    // TITE, MOF, TE.

}
