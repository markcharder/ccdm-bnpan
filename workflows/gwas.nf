process gbz2gfa {

    label 'graphtools'
    label 'process_high'

    storeDir "{$params.gwas_cache}/gbz2gfa"

    input:
    path(gbz)

    output:
    path("graph.gfa"), emit: gfa

    script:
    """
    vg \
        view \
            --threads ${task.cpus} \
            --verbose \
            ${gbz} > graph.gfa
    """

    stub:
    """
    touch graph.gfa
    """

}

process gfa2bin_graph {

    label 'gwastools'
    label 'process_high'

    storeDir "{$params.gwas_cache}/gfa2bin_graph"

    input:
    path(gfa)

    output:
    path("gfa2bin-graph*"), emit: bin

    script:
    """
    gfa2bin \
        graph \
            --gfa ${gfa} \
            --feature node \
            --output gfa2bin-graph \
            --threads ${task.cpus}
    """

    stub:
    """
    touch gfa2bin-graph.{bed,bim,fam}
    """
}

process gfa2bin_cov {

    label 'gwastools'
    label 'process_high'

    storeDir "{$params.gwas_cache}/gfa2bin_cov"

    input:
    tuple path(gfa), path(cov)

    output:
    path("gfa2bin-cov*"), emit: bin

    script:
    """
    gfa2bin \
        graph \
            --gfa ${gfa} \
            --feature node \
            --output graph-${task.index} \
            --threads ${task.cpus}
    """

    stub:
    """
    touch graph-${task.index}.{bed,bim,fam}
    """
}

workflow PACK_CALL {

    take:
    gfa

    main:
    channel
        .fromPath(params.gwas_packs)
        .set { packs }
    
    
    

}
workflow NODE_CALL {

    channel
        .fromPath(params.gwas_graph)
        .set { graph }
    
    if (params.gwas_graphtype == "gbz") {
        gbz2gfa(graph)
            .set { gfa }

    } else if (params.gwas_graphtype == "gfa") {
        gfa = graph

    } else {
        error "Unsupported graph type: ${params.gwas_graphtype}"

    }

    gfa2bin(gfa)

}
