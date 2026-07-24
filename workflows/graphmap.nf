#!/usr/bin/env nextflow

process mash_dist {

    label 'graphmap'
    label 'process_low'
    tag "${variety}"

    input:
    tuple val(variety), path(run1_reads1), path(run1_reads2), path(run2_reads1), path(run2_reads2)

    output:
    tuple val(variety), path("${variety}_combined_R1.fastq.gz"), path("${variety}_combined_R2.fastq.gz"), emit: combined_reads

    script:
    """
    mash sketch -r -m 2 ${run1_reads1}
    mash sketch -r -m 2 ${run1_reads2}
    mash dist ${run1_reads1}.msh ${run2_reads1}.msh
    mash sketch -r -m 2 ${run2_reads1}
    mash sketch -r -m 2 ${run2_reads2}
    mash dist ${run2_reads1}.msh ${run1_reads1}.msh
    """

    stub:
    """
    touch ${variety}_combined_R1.fastq.gz
    touch ${variety}_combined_R2.fastq.gz
    """

}

process kmc_kff {

    label 'graphmap'
    label 'process_high'
    tag "${variety}"
    storeDir "${params.graphmap_cache}/kmc_kff"

    input:
    tuple val(variety), path(reads1), path(reads2)

    output:
    tuple val(variety), path("${variety}.kff"), emit: kff

    script:
    """
    mkdir -p ${variety}_kmc_tmp
    echo -e "${reads1}\\n${reads2}" > ${variety}_kff_input.txt
    kmc \
        -k29 \
        -m${task.memory.toGiga().toInteger()} \
        -t${task.cpus} \
        -okff \
        @${variety}_kff_input.txt \
        ${variety} \
        ${variety}_kmc_tmp 1> ${variety}_kmc.out 2> ${variety}_kmc.err
    """

    stub:
    """
    touch ${variety}.kff
    """

}

process vg_giraffe {

    label 'graphmap'
    label 'process_high'
    tag "${variety}"
    storeDir "${params.graphmap_cache}/vg_giraffe"
    publishDir "${params.graphmap_out_dir}/vg_giraffe", mode: 'copy'

    input:
    tuple val(variety), path(reads1), path(reads2), path(gbz), path(haplotypes), path(kff)

    output:
    path("${variety}_giraffe.gam"), emit: gam

    script:
    """
    vg giraffe \
        --gbz-name ${gbz} \
        --progress \
        --threads ${task.cpus} \
        --fastq-in ${reads1} \
        --fastq-in ${reads2} \
        --haplotype-name ${haplotypes} \
        --kff-name ${kff} \
        --index-basename ccdm-bnpan \
        --sample ${variety} 1> ${variety}_giraffe.gam 2> ${variety}_giraffe.log
    """
    stub:
    """
    touch ${variety}_giraffe.gam
    """

}

process vg_autoindex {

    label 'graphmap'
    label 'autoindex'
    storeDir "${params.graphmap_cache}/vg_autoindex"

    input:
    path(gbz)

    output:
    tuple path("${gbz.baseName}.dist"),
          path("${gbz.baseName}.longread.withzip.min"),
          path("${gbz.baseName}.longread.zipcodes"),
          emit: vg_index_files

    script:
    """
    vg \
        autoindex \
            --threads ${task.cpus} \
            --workflow lr-giraffe \
            --prefix ${gbz.baseName} \
            --gbz ${gbz} \
            --target-mem ${task.memory.toGiga().toInteger()} \
            --verbosity 2
    """

    stub:
    """
    touch ${gbz.baseName}.dist
    touch ${gbz.baseName}.longread.withzip.min
    touch ${gbz.baseName}.longread.zipcodes
    """

}

process vg_giraffe_lr {

    label 'graphmap'
    label 'process_high'
    tag "${variety}"
    publishDir "${params.graphmap_out_dir}/vg_giraffe_lr", mode: 'copy'
    storeDir "${params.graphmap_cache}/vg_giraffe_lr"

    input:
    tuple val(variety), path(reads), path(type), path(dist), path(min), path(zipcodes), path(gbz)

    output:
    path("${variety}_giraffe_lr.gam"), emit: gam_lr

    script:
    """
    vg \
        giraffe \
            --parameter-preset ${type} \
            --gbz-name ${gbz} \
            --fastq-in ${reads} \
            --progress \
            --threads ${task.cpus} \
            --zipcode-name ${zipcodes} \
            --minimizer-name ${min} > ${variety}-lr.gam
    
    vg \
        pack \
            --xg ${gbz} \
            --packs-out ${variety}-lr.pack \
            --gam ${variety}-lr.gam \
            --threads ${task.cpus} \
            --min-mapq 20 > ${variety}-lr.pack
    """

    stub:
    """
    touch ${variety}-lr.gam
    touch ${variety}-lr.pack
    """

}

workflow GRAPHMAP {

/* 
// This is to group reads by variety if needed.
variety_reads_ch
    .groupTuple(by: 0)
    .map { variety, reads1_list, reads2_list ->
        tuple(variety, reads1_list, reads2_list)
    }
    */

    // If we are mapping short reads.
    if (params.graphmap_type == "short") {

        // Read the variety IDs file.
        channel
            .fromPath(params.graphmap_variety_ids)
            .splitCsv(sep: '\t', header: false)
            .map {
                variety, r1, r2 ->
                log.info "${variety} ${r1} ${r2}"
                def prefix = r1.split("\\.")[0]
                tuple("${variety}.${prefix}", file("${params.graphmap_reads_illumina}/${r1}"), file("${params.graphmap_reads_illumina}/${r2}"))
            }
            .set { variety_reads_ch }

            // Generate k-mer frequency files.
            kmc_kff(variety_reads_ch)

            // Combine them with the reads channel and graph+haplotypes.
            channel
                .fromPath(params.graphmap_graph)
                .set { gbz }

            channel
                .fromPath(params.graphmap_haplotypes)
                .set { haplotypes }

            variety_reads_ch
                .combine(gbz)
                .combine(haplotypes)
                .join(kmc_kff.out.kff)
                .set { vg_giraffe_input_ch }

            // Map reads to the graph using vg giraffe.            
            vg_giraffe(vg_giraffe_input_ch)

        }

    else if (params.graphmap_type == "long") {

        channel
            .fromPath(params.graphmap_variety_ids_long)
            .splitCsv(sep: '\t', header: false)
            .map {
                variety, reads, type ->
                tuple(variety, file("${params.graphmap_reads_long}/${reads}"), type)
            }
            .set { variety_reads_ch }
        
        channel
            .fromPath(params.graphmap_graph)
            .set { gbz }
        
        vg_autoindex(gbz)

        variety_reads_ch
            .combine(vg_autoindex.out.vg_index_files)
            .combine(gbz)
            .set { vg_giraffe_lr_input_ch }

        vg_giraffe_lr(vg_giraffe_lr_input_ch)

    }

    else {
        error "Invalid graphmap_type specified in config: ${params.graphmap_type}. Must be 'short' or 'long'."
    }

}
