process minimap_to_ref {

    label 'high'
    label 'minimap2'
    tag "${variety}"

    storeDir "${params.gquals_cache}/minimap2"

    input:
    tuple val(variety), path(genome_fasta), path(reference_fasta)

    output:
    tuple val(variety), path("${variety}-minimap2.paf"), emit: minimap2_out

    script:
    """
    minimap2 -x asm5 ${reference_fasta} ${genome_fasta} > ${variety}-minimap2.paf
    """

    stub:
    """
    touch ${variety}-minimap2.paf
    """

}

process synteny_plot {

    label 'rscript'
    tag "${variety}"

    publishDir "${params.gquals_out_dir}/synteny", mode: 'copy'

    input:
    tuple val(variety), path(minimap2_out)

    output:
    tuple val(variety), path("${variety}-synteny.svg"), emit: synteny_plot_out

    script:
    """
    synteny_dotplot.R ${variety}-minimap2.paf ${variety}
    """

    stub:
    """
    touch ${variety}-synteny.svg
    """
}

process busco_download {

    storeDir "${params.gquals_cache}/busco_lineage"

    output:
    path "${params.gquals_busco_lineage}", emit: busco_lineage_dir

    script:
    """
    wget -O ${params.gquals_busco_lineage}.tar.gz ${params.gquals_busco_lineage_url}
    tar -xvf ${params.gquals_busco_lineage}.tar.gz
    """
    stub:
    """
    mkdir "${params.gquals_busco_lineage}"
    """

}

process busco_aa {

    label 'busco'
    tag "${variety}"
    storeDir "${params.gquals_cache}/busco_aa"
    publishDir "${params.gquals_out_dir}/busco_aa", mode: 'copy'

    input:
    tuple val(variety), path(proteins), path(busco_lineage_dir)

    output:
    tuple val(variety), path("${variety}_busco_aa"), emit: busco

    script:
    """
    busco \
        -i ${proteins} \
        --out ${variety}_busco \
        --mode "protein" \
        --lineage_dataset ${busco_lineage_dir} \
        --cpu ${task.cpus}
    """

    stub:
    """
    mkdir ${variety}_busco
    """

}

workflow GQUALS {

    if (params.gquals_run_busco) {

        busco_download()

        channel
            .fromPath(params.gquals_proteins)
            .map {
                file ->
                def variety = file.baseName.split("_")[0]
                tuple(variety, file)
            }
            .combine(busco_download.out.busco_lineage_dir)
            .set { busco_aa_in_ch }

            busco_aa(busco_aa_in_ch)

    }

    if (params.gquals_run_synteny) {

        channel
            .fromPath(params.gquals_reference)
            .set { reference_ch }
    
        channel
            .fromPath(params.gquals_scaffolded)
            .map {
                file ->
                def variety = file.baseName.split("_")[0]
                tuple(variety, file)
            }
            .set { scaffolded_ch }
    
        scaffolded_ch
            .combine(reference_ch)
            .set { minimap2_input_ch }
    
        minimap_to_ref(minimap2_input_ch)
        synteny_plot(minimap_to_ref.out.minimap2_out)

    }

}