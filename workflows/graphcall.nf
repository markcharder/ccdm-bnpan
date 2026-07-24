// Step 1: calculate snarls.
process vg_snarls {

    label 'process_high'
    label 'graphcall'
    storeDir "${params.graphcall_cache}/vg_snarls"

    input:
    tuple path(gbz)

    output:
    path("${gbz.baseName}.snarls"), emit: snarls

    script:
    """
    vg \
        snarls \
            --threads ${task.cpus} \
            ${gbz} > ${gbz.baseName}.snarls
    """

    stub:
    """
    touch ${gbz.baseName}.snarls
    """

}

// Step 2: filter GAMs.
process vg_filter {

    label 'process_high'
    label 'graphcall'
    tag "${variety}"
    storeDir "${params.graphcall_cache}/vg_filter"

    input:
    tuple val(variety), path(gam), path(gbz)

    output:
    tuple val(variety), path("${variety}-filtered.gam"), emit: gam

    script:
    """
    vg \
        filter \
            --proper-pairs \
            --min-primary 0.9 \
            --frac-score \
            --substitutions \
            --min-end-matches 1 \
            --min-mapq 30 \
            --defray-ends 999 \
            --xg-name ${gbz} \
            --threads ${task.cpus} \
            ${gam} > ${variety}-filtered.gam
    """

    stub:
    """
    touch ${variety}-filtered.gam
    """

}

process vg_pack {

    label 'process_high'
    label 'graphcall'
    tag "${variety}"
    storeDir "${params.graphcall_cache}/vg_pack"

    input:
    tuple val(variety), path(gam), path(gbz)

    output:
    tuple val(variety), path("${variety}.pack"), emit: pack

    script:
    """
    vg \
        pack \
            --xg ${gbz} \
            --gam ${gam} \
            --packs-out ${variety}.pack \
            --threads ${task.cpus}
    """

    stub:
    """
    touch ${variety}.pack
    """

}

process vg_call {

    label 'process_high'
    label 'graphcall'
    tag "${variety}"
    publishDir "${params.graphcall_out_dir}/vg_call", mode: 'copy'

    input:
    tuple val(variety), path(gam), path(pack), path(gbz), path(snarls)

    output:
    tuple val(variety), path("${variety}.vcf"), emit: vcf

    script:
    """
    vg \
        call \
            --genotype-snarls \
            --all-snarls \
            --snarls ${snarls} \
            --sample ${variety} \
            --pack ${pack} \
            --ref-sample ${params.graphcall_vg_refsample} \
            --progress \
            --threads ${task.cpus} \
            --gbz \
            ${gbz} > ${variety}.vcf
    """

    stub:
    """
    touch ${variety}.vcf
    """

}

workflow GRAPHCALL_GAM {

    // Get gbz graph.
    channel
        .fromPath(params.graphcall_graph)
        .set { gbz }

    // Filter snarls.
    vg_snarls(gbz)

    // Get all the input GAMs.
    channel
        .fromPath(params.graphcall_gams)
        .map {
            file ->
                def variety = file.baseName.split('\\.')[0]
                tuple(variety, file)
        }
        .combine(gbz)
        .set { vg_filter_in_ch }

    // Filter the GAMs.
    vg_filter(vg_filter_in_ch)

    // Generate read coverage pack files.
    vg_pack(vg_filter.out.gam)

    // Combine filtered gams with packs, gbz and snarls.
    vg_filter.out.gam
        .join(vg_pack.out.pack)
        .combine(gbz)
        .combine(vg_snarls.out.snarls)
        .set { vg_call_in_ch }

    // Call variants.
    vg_call(vg_call_in_ch)

}

// Hierarchical merging of VCFs using bcftools.
// Step 1: index the VCFs.
process bcftools_index {

    label 'process_medium'
    label 'bcftools'
    storeDir "${params.graphcall_cache}/bcftools_index"

    input:
    path(vcfz) // Takes a single vcf.gz file.

    output:
    tuple path("${vcfz}"), path("${vcfz}.csi"), emit: bcfi // Emits a tuple of vcf.gz and its index.

    script:
    """
    bcftools \
        index \
            --threads ${task.cpus} \
            ${vcfz}
    """

    stub:
    """
    touch ${vcfz}.csi
    """

}

// Step 2: merge chunks of vcfz files.
process bcftools_merge {

    label 'process_medium'
    label 'bcftools'
    storeDir "${params.graphcall_cache}/bcftools_merge"

    input:
    tuple path(vcfzs), path(bcfis)

    output:
    tuple path("bcftools-merged_${task.index}.vcf.gz"), path("bcftools-merged_${task.index}.vcf.gz.csi"), emit: merged_vcf

    script:
    chunk_id = task.index
    println "Got inputs vcfzs: ${vcfzs} bcfis: ${bcfis}"
    println "Chunk id is ${chunk_id}"
    """

    bcftools \
        merge \
            --threads ${task.cpus} \
            --output-type z \
            --output bcftools-merged_${chunk_id}.vcf.gz \
            ${vcfzs.join(' ')}

    bcftools index \
        --threads ${task.cpus} \
        bcftools-merged_${chunk_id}.vcf.gz

    """

    stub:
    """
    touch bcftools-merged_${chunk_id}.vcf.gz
    touch bcftools-merged_${chunk_id}.vcf.gz.csi
    """

}

process bcftools_merge_final {
    
    label 'process_high'
    label 'bcftools'
    publishDir "${params.graphcall_out_dir}/bcftools_merge"

    input:
    tuple path(vcfzs), path(bcfis), val(size) // Size is the number of vcfzs to merge.

    output:
    tuple path("bcftools-merged_f.vcf.gz"), path("bcftools-merged_f.vcf.gz.csi"), emit: merged_vcf

    script:
    println "Got inputs size: ${size} vcfzs: ${vcfzs} bcfis: ${bcfis}"
    if (size > 1) {

        // If there are multiple vcfzs remaining, merge them.
        """

        bcftools \
            merge \
                --threads ${task.cpus} \
                --output-type z \
                --output bcftools-merged_f.vcf.gz \
                ${vcfzs.join(' ')}
        
        bcftools index \
            --threads ${task.cpus} \
            bcftools-merged_f.vcf.gz

        """

    } else {

    // Otherwise rename the vcfz and index to final output name.
        """

        cp \
            ${vcfzs} \
            bcftools-merged_f.vcf.gz

        cp \
            ${bcfis} \
            bcftools-merged_f.vcf.gz.csi

        """

    }

    stub:
    """
    touch bcftools-merged_f.vcf.gz
    touch bcftools-merged_f.vcf.gz.csi
    """

}

workflow GRAPHCALL_BCFMERGE {
    
    // Get a list of vcf.gz files to merge.
    channel
        .fromPath(params.graphcall_bcftools_vcfz)
        .set { vcfz_files }

    // Chunk size provided on command line.
    chunk_size = params.graphcall_bcftools_csize

    // Index the vcf.gz files.
    bcftools_index(vcfz_files)

    // Take output tuple(vcfz, bcfi) and collate into chunks,
    // of tuple(vcfz list, bcfi list).
    bcftools_index.out.bcfi
        .collate(chunk_size)
        .map {
            chunks ->
                def vcfzs = []
                def bcfis = []
                chunks.each { chunk ->
                    vcfzs << chunk[0]
                    bcfis << chunk[1]
                }
            tuple(vcfzs, bcfis)
        }
        .set { merge1_ch }

    // First, merge chunks. Returns a tuple of 2 lists if there were multiple chunks.
    merge1 = bcftools_merge(merge1_ch)

    // Take the tuple of 2 lists and calculate its size, then pass to final merge.
    merge1.merged_vcf
        .collect(flat: false)
        .map { tuples ->
            def vcfzs = []
            def bcfis = []
            tuples.each { t ->
                vcfzs << t[0]
                bcfis << t[1]
            }
            def size = vcfzs.size()
            tuple(vcfzs, bcfis, size)
        }
        .set { merge2_ch }

    // Run final merge.
    bcftools_merge_final(merge2_ch)

}
