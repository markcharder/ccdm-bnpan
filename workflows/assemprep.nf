#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ragtag_scaffold {

    label 'process_medium'
    label 'medium'
    label 'assemprep'

    tag "${sample}"

    storeDir "${params.assemprep_cache}/${sample}_ragtag"

    input:
    tuple val(sample), path(query)
    path reference

    output:
    tuple val(sample), path("${sample}_ragtag/"), emit: ragtag_out

    script:
    """
    ragtag.py scaffold \\
        -t ${task.cpus} \\
        -o ${sample}_ragtag \\
        -u \\
        ${reference} \\
        ${query}
    """

    stub:
    """
    mkdir -p ${sample}_ragtag
    touch ${sample}_ragtag/ragtag.scaffold.fasta
    """

}

process split_scaffolds_by_genome {

    label 'small'
    label 'assemprep'

    tag "${sample}"

    publishDir "${params.assemprep_out_dir}/${sample}_split", mode: 'copy'

    input:
    tuple val(sample), path(ragtag_out)
    path genome_tsv

    output:
    tuple val(sample), path("${sample}_A.fasta"),        emit: a_genome
    tuple val(sample), path("${sample}_C.fasta"),        emit: c_genome
    tuple val(sample), path("${sample}_unplaced.fasta"), emit: unplaced

    script:
    """
    split_scaffolds_by_genome.py \\
        --fasta ${ragtag_out}/ragtag.scaffold.fasta \\
        --tsv   ${genome_tsv} \\
        --sample ${sample}
    """

    stub:
    """
    touch ${sample}_A.fasta ${sample}_C.fasta ${sample}_unplaced.fasta
    """

}

workflow ASSEMPREP {

    reference = file(params.assemprep_reference)
    genome_tsv = file(params.assemprep_tsv)

    channel
        .fromPath("${params.assemprep_assemblies_dir}/*.{fa,fa.gz,fasta,fasta.gz,fna,fna.gz}")
        .map { f -> [ f.simpleName, f ] }
        .filter { _sample, f -> f.toRealPath() != reference.toRealPath() }
        .set { genomes_ch }

    ragtag_scaffold(genomes_ch, reference)

    split_scaffolds_by_genome(ragtag_scaffold.out.ragtag_out, genome_tsv)

}

