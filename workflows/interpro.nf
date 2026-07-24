#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process download_interpro_db {

    label 'small' // This is a small, serial process that requires minimal memory. Don't waste the HPC.

    storeDir params.interpro_db_cache

    output:
    path "${params.interpro_db_name}"

    script:
    """
    wget -O "${params.interpro_db_name}.tar.gz" "${params.interpro_db_url}"
    tar -xvf "${params.interpro_db_name}.tar.gz"
    """

    stub:
    """
    mkdir -p "${params.interpro_db_name}"
    """

}

process split_fasta {

    label 'seqkit'
    tag "${variety}"

    input:
    tuple val(variety), path(protein_fasta)

    output:
    tuple val(variety), path("${protein_fasta}.split/*"), emit: split_fasta

    script:
    def is_gzipped = protein_fasta.name.endsWith('.gz')
    def input_file = is_gzipped ? "braker.aa" : protein_fasta
    def decompress = is_gzipped ? "gunzip -c ${protein_fasta} > braker.aa" : ""

    """
    ${decompress}
    seqkit split2 --by-size "${params.interpro_chunk_size}" "${input_file}"
    """

}

process interproscan {

    label 'interproscan'
    tag "${variety}"

    input:
    tuple val(variety), path(protein_chunk), val(db_ready)

    output:
    tuple val(variety), path("${variety}-${protein_chunk}_interpro.tsv"), emit: interpro_results

    script:
    """
    if [ ! -d "/opt/interproscan/data" ]; then
        echo "I can't find the interproscan database directory from inside the interproscan container..."
        echo "Process 'interproscan' says goodbye, my friend."
        exit 1
    fi
    /opt/interproscan/interproscan.sh \
        --input ${protein_chunk} \
        --formats tsv \
        -o ${variety}-${protein_chunk}_interpro.tsv \
        --cpu ${task.cpus} \
        --outfile ${variety}-${protein_chunk}_interpro.tsv \
        --tempdir interpro_temp/ \
        --goterms \
        --pathways \
        --seqtype p \
        --verbose
    """

}

process combine_interpro {

    label 'small'
    tag "${variety}"
    publishDir "${params.interpro_out_dir}", mode: 'copy'

    input:
    tuple val(variety), path(interpro_files)

    output:
    path "${variety}_interpro_combined.tsv"

    script:
    """
    cat ${interpro_files.join(' ')} > ${variety}_interpro_combined.tsv
    """

}

workflow INTERPRO {

    take:
    braker_proteins_ch

    main:
    download_interpro_db()
    split_fasta(braker_proteins_ch)
    split_fasta.out.split_fasta
        .transpose()
        .combine(download_interpro_db.out)
        .set { protein_chunks_ch }
    interproscan(protein_chunks_ch)
    interproscan.out.interpro_results
        .groupTuple(by: 0)
        .map { variety, interpro_files ->
            tuple(variety, interpro_files)
        }
        .set { interpro_grouped_ch }
    combine_interpro(interpro_grouped_ch)
} 