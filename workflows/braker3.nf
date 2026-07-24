#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def braker3HelpMessage() {

    log.info"""

    BRAKER3 WORKFLOW
    ============================
    This workflow performs genome annotation using BRAKER3. It runs annotation on one genome first
    and uses the parameters from this run to annotate the remaining genomes.

    Optional parameters:


        GENERAL:

        --help                    Show this help message and exit.
        
        --braker3_masked_genomes  Directory containing the masked genomes to annotate.
                                  Default: data/processed/rnaseq/masked_genomes_full/

        --braker3_bam_dir         Directory containing the RNA-Seq alignments in BAM format.
                                  Default: data/processed/rnaseq/alignments/

        --braker3_out_dir         Output for braker3 annotations.
                                  Default : data/processed/braker3/
        
        --braker3_prot_db_url     URL for downloading odb for braker3 to use as the input protein database.
                                  Default: https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Viridiplantae.fa.gz

        --braker3_prot_db_name    Name of the file to save the downloaded protein database to.
                                  Default: Viridiplantae.fa.gz

        --braker3_cache           Directory for caching braker3 runs, including augustus config for the first one.
                                  Default: .braker3/

        --braker3_prot_db_cache   Directory for caching the protein database for braker3.
                                  Default: data/databases/

        --braker3_first           Name of variety for the first braker3 run.
                                  Default: 'Westar'

        --braker3_skip_rest       Whether or not to skip braker3 for genomes beyond the first run.
                                  Default: false

        --braker3_species_name    The species name to use for braker3. 
                                  Default: "bnpan"

""".stripIndent()

}

process download_protein_db {

    output:
    path "${params.braker3_prot_db_name}".replace('.gz',''), emit: protein_db

    script:
    """
    wget -O "${params.braker3_prot_db_name}" "${params.braker3_prot_db_url}"
    gunzip "${params.braker3_prot_db_name}"
    """

    stub:
    """
    touch "${params.braker3_prot_db_name}".replace('.gz','')
    """

}

process braker3_one_genome {

    label 'braker3'
    tag "${variety}"

    publishDir "${params.braker3_out_dir}/${variety}_braker3", mode: 'copy'
    storeDir "${params.braker3_cache}/${variety}_braker3_initial"

    input:
    tuple val(variety),
        val(bams_list),
        path(bams),
        path(masked_genome),
        path(protein_db)

    output:
    tuple val(variety), path("${variety}_braker3/braker.aa"), emit: annot_proteins
    tuple val(variety), path("${variety}_braker3/braker.gtf"), emit: annot_gtf
    tuple val(variety), path("${variety}_braker3/braker.codingseq"), emit: annot_codingseq
    path "augustus_config", emit: augustus_config

    script:
    """

    # Record current directory.
    export WORKDIR=\$(pwd)

    # Create a /tmp directory to run this on tmpfs.
    export TMPDIR=/tmp/braker_${variety}_${task.index}_\$(date +%s)
    mkdir -p \$TMPDIR/

    # Export AUGUSTUS config path to a writable location.
    export AUGUSTUS_CONFIG_PATH=\$WORKDIR/augustus_config
    mkdir -p \$AUGUSTUS_CONFIG_PATH
    
    # Copy existing contents of augustus config to the new location.
    cp -r /opt/Augustus/config/* \$AUGUSTUS_CONFIG_PATH/ 
    
    # Stage data to tmpdir for fast IO.
    cp ${masked_genome} \$TMPDIR/
    cp ${protein_db} \$TMPDIR/
    cp ${bams} \$TMPDIR/

    # Change to tmpdir.
    cd \$TMPDIR/
    
    # Run BRAKER3.
    braker.pl \
        --AUGUSTUS_CONFIG_PATH=\"\$AUGUSTUS_CONFIG_PATH\" \
        --AUGUSTUS_BIN_PATH=/opt/Augustus/bin \
        --AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts \
        --genome ${masked_genome} \
        --prot_seq ${protein_db} \
        --bam ${bams_list} \
        --species ${params.braker3_species_name} \
        --verbosity 4 \
        --workingdir ${variety}_braker3 \
        --threads ${task.cpus} &
    
    BRAKER_PID=\$!
    mkdir -p \$WORKDIR/tmp_errors/

    while kill -0 "\$BRAKER_PID" 2> /dev/null; do
    
        if [ -f "${variety}_braker3/braker.log" ]; then
            cp "${variety}_braker3/braker.log" \
            \$WORKDIR/

        fi

        if ls ${variety}_braker3/errors/* 1> /dev/null 2>&1; then
                cp -r ${variety}_braker3/errors/* \
                \$WORKDIR/tmp_errors/

            fi

        sleep 10

    done

    wait \$BRAKER_PID
    BRAKER_EXIT=\$?

    # Copy results back to workdir
    mv \$TMPDIR/${variety}_braker3/ \$WORKDIR/
    rm -rf \$TMPDIR
    """

    stub:
    """
    mkdir -p "${variety}_braker3"
    mkdir -p "augustus_config/species/${params.braker3_species_name}"
    touch "${variety}_braker3/braker.aa"
    touch "${variety}_braker3/braker.gtf"
    touch "${variety}_braker3/braker.codingseq"
    """

}

process braker3 {

    label 'braker3'
    tag "${variety}"
        
    publishDir "${params.braker3_out_dir}/${variety}_braker3", mode: 'copy'
    storeDir "${params.braker3_cache}/${variety}_braker3"

    input:
    tuple val(variety),
        val(bams_list),
        path(bams),
        path(masked_genome),
        path(protein_db),
        path(augustus_config)
    
    output:
    tuple val(variety), path("${variety}_braker3/braker.aa"), emit: annot_proteins
    tuple val(variety), path("${variety}_braker3/braker.gtf"), emit: annot_gtf
    tuple val(variety), path("${variety}_braker3/braker.codingseq"), emit: annot_codingseq

    script:
    """

    # We get the current working directory of the process, set a space aside in tmp
    # and set a unique AUGUSTUS config path for this task.
    export WORKDIR=\$(pwd)
    export TMPDIR=/tmp/braker_${task.index}_\$(date +%s)
    mkdir -p augconfwd
    cp -r "${augustus_config}"/* augconfwd/
    export AUGUSTUS_CONFIG_PATH=\$WORKDIR/augconfwd

    # We make the tmpdir.
    mkdir -p \$TMPDIR/

    # Stage the data to the tmpdir for fast IO.
    cp ${masked_genome} \$TMPDIR/
    cp ${protein_db} \$TMPDIR/
    cp ${bams} \$TMPDIR/

    # Change ourselves to the tmpdir because it's nice in there.
    cd \$TMPDIR/

    # Run braker.pl in the tmpdir so that prothint has fast IO.
    braker.pl \
        --AUGUSTUS_CONFIG_PATH=\"\$AUGUSTUS_CONFIG_PATH\" \
        --AUGUSTUS_BIN_PATH=/opt/Augustus/bin \
        --AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts \
        --useexisting \
        --genome ${masked_genome} \
        --prot_seq ${protein_db} \
        --bam ${bams_list} \
        --verbosity 4 \
        --workingdir ${variety}_braker3 \
        --species ${params.braker3_species_name} \
        --threads ${task.cpus} &
    
    BRAKER_PID=\$!
    mkdir -p \$WORKDIR/tmp_errors/

    while kill -0 "\$BRAKER_PID" 2> /dev/null; do
    
        if [ -f "${variety}_braker3/braker.log" ]; then
            cp "${variety}_braker3/braker.log" \
            \$WORKDIR/

        fi

        if ls ${variety}_braker3/errors/* 1> /dev/null 2>&1; then
                cp -r ${variety}_braker3/errors/* \
                \$WORKDIR/tmp_errors/

            fi

        sleep 10

    done

    wait \$BRAKER_PID
    BRAKER_EXIT=\$?

    # Copy the results back to the workdir.
    # These should be available in the output channel for the next process.
    mv \$TMPDIR/${variety}_braker3/ \$WORKDIR
    rm -rf \$TMPDIR

    """

    stub:
    """
    mkdir -p "${variety}_braker3"
    touch "${variety}_braker3/braker.gtf"
    touch "${variety}_braker3/braker.aa"
    touch "${variety}_braker3/braker.gff3"
    """

}

process braker3_rename {

    tag "${variety}"
    publishDir "${params.braker3_out_dir}/braker3_renamed", mode: 'copy'

    input:
    tuple val(variety), path(braker3_aa)

    output:
    tuple val(variety), path("${variety}_braker3.fa"), emit: braker3_aa_renamed

    script:
    """
    cp ${braker3_aa} ${variety}_braker3.fa
    """

    stub:
    """
    touch ${variety}_braker3.fa
    """

}

workflow BRAKER3 {

    main:
    if (params.help) {
        braker3HelpMessage()
        System.exit(0)
    }

    // If necessary, start downloading the protein database.
    download_protein_db()
 
    // Get bams for annotation.
    channel
        .fromPath(
            "${params.braker3_bam_dir}/*.bam"
            )
        .map {
            bams ->
            def variety = bams.baseName.split('_')[0]
            tuple(variety, bams)
        }
        .groupTuple() 
        .set { bams_by_variety }

    // Get repeat masked genome for annotation.
    channel
        .fromPath("${params.braker3_masked_genomes}/*.{fa,fasta}")
        .map {
            file ->
            def variety = file.baseName.split('_')[0]
            tuple(variety, file)
        }
        .set { genomes_by_variety }
 
    bams_by_variety
        .join(genomes_by_variety)
        .map {
            variety, bams, genome ->
            def bams_list = bams.join(",")
            tuple(variety, bams_list, bams, genome)
        }
        .combine(download_protein_db.out.protein_db)
        .set { braker_input_ch }

    braker_input_ch
        .filter {
            variety, _bams_list, _bams, _genome, _protein_db -> 
                variety == params.braker3_first
            }
        .set { braker3_first_ch }

    braker3_one_genome(braker3_first_ch)

    if (!params.braker3_skip_rest) {
        braker_input_ch
            .filter {
                variety, _bams_list, _bams, _genome, _protein_db -> 
                    variety != params.braker3_first
                }
            .map {
                variety, bams_list, bams, genome, protein_db ->
                tuple(variety, bams_list, bams, genome, protein_db)
            }
            .combine(braker3_one_genome.out.augustus_config)
            .set { braker3_rest_ch }

        braker3(braker3_rest_ch)
/*
    WORKING ON THIS:
        if (params.braker3_rename_aa) {
            braker3_one_genome.out.annot_proteins 
                    
                }

        }
*/
    }

}