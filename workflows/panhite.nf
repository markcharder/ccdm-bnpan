/*  This was forked from the panHiTE workflow and modified to streamline
    integration.

    Processes are:

        1) pan_preprocess_genomes
        2) pan_run_hite_single
        3) merge_terminal_te
        4) pan_remove_redundancy
        5) pan_recover_low_copy_TEs
        6) pan_merge_TE_recover
        7) pan_annotate_genomes
        8) pan_summarize_tes
        9) pan_gene_te_relation
        10) pan_generate_bam_for_RNA_seq
        11) pan_detect_de_genes
        12) pan_split_genome
        13) pan_annotate_chunk
        14) pan_merge_annotations
        15) pan_recover_split_annotation

*/

// 1 
process pan_preprocess_genomes {

    label 'panhite'
    storeDir "${params.out_dir}/pan_preprocess_genomes"

    input:
    path genome_list
    path genes_dir
    path RNA_dir
    path pan_genomes_dir

    output:
    path "genome_metadata.json"

    script:
    """
    echo "Running preprocessing"
    pan_preprocess_genomes.py --genome_list ${genome_list} --genes_dir ${genes_dir} \
    --RNA_dir ${RNA_dir} --pan_genomes_dir ${pan_genomes_dir} > preprocess_genomes.log 2>&1
    """

    stub:
    """
    touch genome_metadata.json
    """
}

process pan_run_hite_single {

    label 'panhite'
    storeDir "${params.out_dir}/pan_run_hite_single/${genome_name}"

    input:
    tuple val(genome_name), val(raw_name), path(reference), val(te_type), val(miu), val(debug)

    output:
    tuple val(genome_name), path("intact_LTR.list"), emit: ch_intact_ltr_list
    path "confident_ltr.terminal.fa", emit: ch_ltr_terminal
    path "confident_ltr.internal.fa", emit: ch_ltr_internal
    path "confident_ltr_cut.fa", emit: ch_ltr_cut
    path "filtered_single_copy_LTRs.txt", emit: ch_filtered_single_ltr
    path "intact_LTR.fa", emit: ch_intact_ltr
    path "intact_LTR.fa.classified", emit: ch_intact_classified_ltr
    path "confident_helitron.fa", emit: ch_helitron
    path "confident_non_ltr.fa", emit: ch_non_ltr
    path "confident_other.fa", emit: ch_other
    path "confident_tir.fa", emit: ch_tir
    path "confident_TE.cons.fa", emit: ch_te
    path "chr_name.map", emit: chr_name_map
    tuple val(genome_name), path("tir_low_copy.fa"), path("helitron_low_copy.fa"), path("non_ltr_low_copy.fa"), emit: ch_low_copy_files

    script:
    cores = task.cpus
    """
    pan_run_hite_single.py --genome_name ${genome_name} --reference ${reference} --threads ${cores} \
    --te_type ${te_type} --miu ${miu} --work_dir ${params.work_dir} --debug ${debug} --shared_prev_TE ${params.shared_prev_TE} > ${genome_name}.run_hite_single.log 2>&1
    """

    stub:
    """
    # touch ing_story_of_all_expected_output_files

    touch confident_ltr.terminal.fa
    touch confident_ltr.internal.fa
    touch confident_ltr_cut.fa
    touch filtered_single_copy_LTRs.txt
    touch intact_LTR.fa
    touch intact_LTR.fa.classified
    touch confident_helitron.fa
    touch confident_non_ltr.fa
    touch confident_other.fa
    touch confident_tir.fa
    touch confident_TE.cons.fa
    touch chr_name.map
    touch tir_low_copy.fa
    touch helitron_low_copy.fa
    touch non_ltr_low_copy.fa
    """
}

process merge_terminal_te {

    label 'panhite'
    storeDir "${params.out_dir}/merge_terminal_te"

    input:
    path terminal
    path tir
    path helitron
    path other
    path non_ltr

    output:
    path "pan_terminal_te.merge.fa"

    script:
    """
    cat ${terminal} ${tir} ${helitron} ${other} ${non_ltr} > pan_terminal_te.merge.fa
    """

    stub:
    """
    # touch ed_by_the_effort_that_went_into_the_amusing_comments

    touch pan_terminal_te.merge.fa
    """
}

process pan_remove_redundancy {

    label 'panhite'
    storeDir "${params.out_dir}/pan_remove_redundancy"

    input:
    file merge_te_file

    output:
    path "panTE.fa"

    script:
    cores = task.cpus
    """
    pan_remove_redundancy.py --merge_te_file ${merge_te_file} --threads ${cores} > pan_remove_redundancy.log 2>&1
    """

    stub:
    """
    # touch down

    touch panTE.fa
    """
}

process pan_recover_low_copy_TEs{

    label 'panhite'
    storeDir "${params.out_dir}/pan_recover_low_copy_TEs/${genome_name}"

    input:
    tuple val(genome_name), path(tir_low_copy), path(helitron_low_copy), path(non_ltr_low_copy), path(panTE_lib)

    output:
    path "panTE.recover.fa.classified", emit: ch_recover_TEs
    path "get_copies", emit: chr_copies
    path "real_*.fa.cons", emit: chr_raw_real

    script:
    cores = task.cpus
    """
    pan_recover_low_copy_TEs.py --genome_name ${genome_name} --threads ${cores} --tir_low_copy ${tir_low_copy} \
    --helitron_low_copy ${helitron_low_copy} --non_ltr_low_copy ${non_ltr_low_copy} \
    --panTE_lib ${panTE_lib} --genome_list ${params.genome_list} \
    --pan_genomes_dir ${params.pan_genomes_dir} > pan_recover_low_copy_TEs.log 2>&1
    """

    stub:
    """
    # enough

    touch panTE.recover.fa.classified
    touch get_copies
    touch real_really.fa.cons
    """
}

process pan_merge_TE_recover {

    label 'panhite'
    storeDir "${params.out_dir}/pan_merge_TE_recover"

    input:
    path panTE_lib
    path panTE_recover_lib

    output:
    path "panTE.merge_recover.fa", emit: ch_panTE_merge
    path "panTE.merge_recover.fa.clstr", emit: ch_clstr

    script:
    cores = task.cpus
    """
    cat ${panTE_lib} ${panTE_recover_lib} > panTE.merge_recover.redundant.fa
    cd-hit-est -aS 0.95 -aL 0.95 -c 0.8 -d 0 -G 0 -g 1 -A 80 -i panTE.merge_recover.redundant.fa \
    -o panTE.merge_recover.fa -T ${cores} -M 0 > pan_merge_TE_recover.log 2>&1
    """

    stub:
    """
    touch panTE.merge_recover.fa
    touch panTE.merge_recover.fa.clstr
    """
}

process pan_annotate_genomes {

    label 'panhite'
    storeDir "${params.out_dir}/pan_annotate_genomes/${genome_name}"

    input:
    tuple val(genome_name), path(reference), path(panTE_lib)

    output:
    tuple val(genome_name), path("${genome_name}.gff"), path("${genome_name}.tbl"), path("${genome_name}.full_length.gff"), path("${genome_name}.full_length.copies"), emit: annotate_out

    script:
    cores = task.cpus
    """
    pan_annotate_genome.py --threads ${cores} --panTE_lib ${panTE_lib} --reference ${reference} \
    --genome_name ${genome_name} > ${genome_name}.annotate_genomes.log 2>&1
    """
    stub:
    """
    touch ${genome_name}.gff
    touch ${genome_name}.tbl
    touch ${genome_name}.full_length.gff
    touch ${genome_name}.full_length.copies
    """
}

process pan_summarize_tes {

    label 'panhite'
    storeDir "${params.out_dir}/pan_summarize_tes"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.pdf"

    input:
    path genome_info_json
    path pan_genomes_dir
    path panTE_lib
    val softcore_threshold

    output:
    path "TE_summary.pdf", emit: ch_te_summary
    path "panHiTE.CorePan_fitmodel.pdf", emit: ch_corepan_model, optional: true
    path "panHiTE.CorePan_fitsmooth.pdf", emit: ch_corepan_smodel, optional: true
    path "panHiTE_fl.CorePan_fitmodel.pdf", emit: ch_corepan_fl_model, optional: true
    path "panHiTE_fl.CorePan_fitsmooth.pdf", emit: ch_corepan_fl_smodel, optional: true
    path "*.txt", emit: ch_corepan_txt, optional: true
    path "TE_PAV.tsv", emit: ch_pav
    path "full_length_TE_PAV.tsv", emit: ch_full_length_pav
    path "Full_length_TEs_Ratio.json", emit: ch_fl_te_ratio
    path "TEs_Ratio.json", emit: ch_te_ratio
    path "Full_length_TE_Coverage.json", emit: ch_fl_te_coverage
    path "TE_Coverage.json", emit: ch_te_coverage
    path "Full_length_TE_Classes_Ratio.json", emit: ch_fl_te_classes_ratio
    path "TE_Classes_Ratio.json", emit: ch_te_classes_ratio
    path "Full_length_TE_Classes_Coverage.json", emit: ch_fl_te_classes_coverage
    path "TE_Classes_Coverage.json", emit: ch_te_classes_coverage
    path "intact_LTR_insert_time.csv", emit: ch_ltr_insert_time

    script:
    """
    pan_summary_TEs.py --genome_info_json ${genome_info_json} --pan_genomes_dir ${pan_genomes_dir} \
    --panTE_lib ${panTE_lib} --softcore_threshold ${softcore_threshold} > summarize_tes.log 2>&1
    """

    stub:
    """
    touch TE_summary.pdf
    touch panHiTE.CorePan_fitmodel.pdf
    touch panHiTE.CorePan_fitsmooth.pdf
    touch panHiTE_fl.CorePan_fitmodel.pdf
    touch panHiTE_fl.CorePan_fitsmooth.pdf
    touch txt.txt
    touch TE_PAV.tsv
    touch full_length_TE_PAV.tsv
    touch full_length_TEs_Ratio.json
    touch TEs_Ratio.json
    touch Full_length_TE_Coverage.json
    touch TE_Coverage.json
    touch Full_length_TE_Classes_Ratio.json
    tuoch TE_Classes_Ratio.json
    touch Full_length_TE_Classes_Coverage.json
    touch TE_Classes_Coverage.json
    touch intact_LTR_insert_time.csv
    """
}


process pan_gene_te_relation {

    label 'panhite'
    storeDir "${params.out_dir}/pan_gene_te_relation"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "gene_te_associations.tsv"

    input:
    path genome_info_json

    output:
    path "gene_te_associations.tsv"

    script:
    """
    pan_gene_te_relation.py --genome_info_json ${genome_info_json} > pan_gene_te_relation.log 2>&1
    """
    stub:
    """
    touch gene_te_associations.tsv
    """
}


process pan_generate_bam_for_RNA_seq {

    label 'panhite'
    storeDir "${params.out_dir}/pan_generate_bam_for_RNA_seq"

    input:
    tuple val(genome_name), path(reference), val(RNA_dir), val(RNA_seq)

    output:
    tuple val(genome_name), path("${genome_name}.output.sorted.bam"), emit: bam_out

    script:
    cores = task.cpus
    """
    pan_generate_bam_for_RNA_seq.py --genome_name ${genome_name} --reference ${reference} \
    --RNA_seq '${RNA_seq}' --RNA_dir ${RNA_dir} --threads ${cores}  > ${genome_name}.pan_generate_bam_for_RNA_seq.log 2>&1

    sync
    """
    stub:
    """
    touch "${genome_name}.output.sorted.bam"
    """
}

process pan_detect_de_genes {

    label 'panhite'
    storeDir "${params.out_dir}/pan_detect_de_genes"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.tsv"

    input:
    path genome_info_for_bam_json
    path gene_te_associations
    val RNA_dir

    output:
    path "DE_genes_from_TEs.tsv", emit: ch_de_genes, optional: true
    path "all_gene_TEs_details.tsv", emit: ch_all_genes, optional: true
    path "DE_genes_from_TEs.pdf", emit: ch_de_genes_pdf, optional: true
    path "TE_express.table", emit: ch_te_express, optional: true
    path "gene_express.table", emit: ch_gene_express, optional: true

    script:
    cores = task.cpus
    """
    pan_detect_de_genes.py --genome_info_for_bam_json ${genome_info_for_bam_json} --gene_te_associations ${gene_te_associations} \
    --RNA_dir ${RNA_dir} --threads ${cores} > pan_detect_de_genes.log 2>&1
    """
    stub:
    """
    touch DE_genes_from_TEs.tsv
    touch all_gene_TEs_details.tsv
    touch DE_genes_from_TEs.pdf
    touch TE_express.table
    touch gene_express.table
    """
}

process pan_split_genome {

    label 'panhite'
    storeDir "${params.out_dir}/pan_split_genome/${genome_name}"

    input:
    tuple val(genome_name), path(reference)

    output:
    tuple val(genome_name), path("genome.cut*.fa"), emit: ch_chunks

    script:
    """
    split_genome_chunks.py -g ${reference} --chrom_seg_length ${params.chrom_seg_length} \
    --chunk_size ${params.chunk_size} > ${genome_name}.pan_split_genome.log 2>&1
    """

    stub:
    """
    touch genome.cut1.fa
    touch genome.cut2.fa
    """
}

process pan_annotate_chunk {

    label 'panhite'
    storeDir "${params.out_dir}/pan_annotate_chunk/${genome_name}"

    input:
    tuple val(genome_name), path(chunk_fasta), path(panTE_lib)

    output:
    tuple val(genome_name), path("${genome_name}_${ref_index}.gff"), path("${genome_name}_${ref_index}.full_length.gff"), path("${genome_name}_${ref_index}.full_length.copies"), emit: ch_annotated_chunks

    script:
    (_full, ref_index) = (chunk_fasta =~ /genome.cut(\d+)\.fa/)[0]
    """
    pan_annotate_genome.py --threads ${task.cpus} --panTE_lib ${panTE_lib} \
    --reference ${chunk_fasta} --genome_name ${genome_name}_${ref_index} > ${genome_name}_${ref_index}.pan_annotate_chunk.log 2>&1
    """
    stub:
    """
    touch ${genome_name}_${ref_index}.gff"
    touch ${genome_name}_${ref_index}.full_length.gff"
    touch ${genome_name}_${ref_index}.full_length.copies"
    """
}

process pan_merge_annotations {

    label 'panhite'
    storeDir "${params.out_dir}/pan_merge_annotations/${genome_name}"

    input:
    tuple val(genome_name), path(annotated_chunks), path(full_length_annotated_chunks), path(full_length_copies_chunks)

    output:
    tuple val(genome_name), path("${genome_name}_merged.gff"), path("${genome_name}_merged.full_length.gff"), path("${genome_name}_merged.full_length.copies"), emit: ch_merged_annotations

    script:
    """
    cat ${annotated_chunks} > ${genome_name}_merged.gff
    cat ${full_length_annotated_chunks} > ${genome_name}_merged.full_length.gff
    awk '1; END {print ""}' ${full_length_copies_chunks} > ${genome_name}_merged.full_length.copies
    """
    stub:
    """
    touch ${genome_name}_merged.gff
    touch ${genome_name}_merged.full_length.gff
    touch ${genome_name}_merged.full_length.copies
    """
}

process pan_recover_split_annotation {

    label 'panhite'
    storeDir "${params.out_dir}/pan_annotate_genomes/${genome_name}"

    input:
    tuple val(genome_name), path(input_gff), path(full_length_input_gff), path(full_length_copies), path(genome_path), path(te_lib)

    output:
    tuple val(genome_name), path("${genome_name}.sorted.gff"), path("${genome_name}.sorted.gff.tbl"), path("${genome_name}.full_length.sorted.gff"), path("${genome_name}_merged.full_length.copies"), emit: annotate_out

    script:
    """
    recover_split_annotation.py ${input_gff} ${genome_name}.gff
    bedtools sort -i ${genome_name}.gff > ${genome_name}.sorted.gff
    get_summary_count.sh ${te_lib} ${genome_name}.sorted.gff ${genome_path}

    recover_split_annotation.py ${full_length_input_gff} ${genome_name}.full_length.gff
    bedtools sort -i ${genome_name}.full_length.gff > ${genome_name}.full_length.sorted.gff
    get_summary_count.sh ${te_lib} ${genome_name}.full_length.sorted.gff ${genome_path} -gff3
    """
    stub:
    """
    # This is the last stub, aaaaaaaaaaah yes!
    # touch ed_by_the_effort_put_into_making_stubs
    touch ${genome_name}.sorted.gff
    touch ${genome_name}.sorted.gff.tbl
    touch ${genome_name}.full_length.sorted.gff
    touch ${genome_name}_merged.full_length.copies
    """
}

workflow panHiTE_workflow {

    main:

    // Exit if the wrong TE type is given.
    def all_te_types = ['ltr', 'tir', 'helitron', 'non-ltr', 'all']

    if (!all_te_types.contains(params.te_type)) {
        error "Invalid TE type: ${params.te_type}. Please choose from ${all_te_types}"
    }

    genome_metadata_out = pan_preprocess_genomes(
        params.genome_list,
        params.genes_dir,
        params.RNA_dir,
        params.pan_genomes_dir
        )

    genome_info_list = genome_metadata_out
        .flatMap { json_file ->
            def genome_data = file(json_file).text
            def genome_json = new groovy.json.JsonSlurper().parseText(genome_data)
            genome_json.genome_info.collect { it -> [it.genome_name, it.raw_name, it.reference, it.gene_gtf, it.RNA_seq] }
        }

    genome_info_list.map { genome_name, raw_name, reference, _gene_gtf, _RNA_seq ->
            [genome_name, raw_name, file(reference), params.te_type, params.miu, params.debug]
        }.set { hite_input_channel }

    hite_out = pan_run_hite_single(hite_input_channel)
    all_te = hite_out.ch_te.collectFile(name: "${params.out_dir}/pan_te.tmp.fa")
    intact_ltr_list_channel = hite_out.ch_intact_ltr_list
    low_copy_files_channel = hite_out.ch_low_copy_files

    panTE_lib = pan_remove_redundancy(all_te)
    panTE_lib = panTE_lib.collectFile(name: "${params.out_dir}/panTE.fa")

    _recover_input = low_copy_files_channel.combine(panTE_lib).set { recover_input_channel }

    recover_out = pan_recover_low_copy_TEs(recover_input_channel)
    panTE_recover_lib = recover_out.ch_recover_TEs
    panTE_recover_lib = panTE_recover_lib.collectFile(name: "${params.out_dir}/panTE.recover.fa")

    merge_out = pan_merge_TE_recover(panTE_lib, panTE_recover_lib)
    panTE_merge_lib = merge_out.ch_panTE_merge
    panTE_merge_lib = panTE_merge_lib.collectFile(name: "${params.out_dir}/panTE.merge_recover.fa")
    genome_info_list.map { genome_name, _raw_name, reference, _gene_gtf, _RNA_seq ->
            [genome_name, reference]
        }.set { split_genome_input_channel }
    split_genome_out = pan_split_genome(split_genome_input_channel)

    split_genome_out.ch_chunks.flatMap { genome_name, chunk_files ->
        chunk_files = chunk_files instanceof List ? chunk_files : [chunk_files]
        chunk_files.collect { chunk_file ->
            [genome_name, chunk_file]
            }
        }
        .combine(panTE_merge_lib)
        .set { annotate_input_channel }
    annotate_chunk_out = pan_annotate_chunk(annotate_input_channel)
    annotate_chunk_out.ch_annotated_chunks
        .groupTuple(by: 0)
        .set { grouped_chunks }
    merge_annotations_out = pan_merge_annotations(grouped_chunks)
    merge_annotations_out.ch_merged_annotations.join(split_genome_input_channel).combine(panTE_merge_lib).set { recover_split_annotation_input_channel }
    annotate_out = pan_recover_split_annotation(recover_split_annotation_input_channel)
    genome_info_list.join(annotate_out).map { genome_info ->
        [
            genome_info[0],
            genome_info[1],
            genome_info[2],
            genome_info[3],
            genome_info[4],
            genome_info[5].toString(),
            genome_info[6].toString(),
            genome_info[7].toString(),
            genome_info[8].toString()
        ]
    }.join(intact_ltr_list_channel).map { genome_info ->
        [
            genome_name: genome_info[0],
            raw_name   : genome_info[1],
            reference  : genome_info[2],
            gene_gtf   : genome_info[3],
            RNA_seq    : genome_info[4],
            TE_gff        : genome_info[5],
            full_length_TE_gff: genome_info[7],
            full_length_copies: genome_info[8],
            intact_LTR_list: genome_info[9].toString()
        ]
    }.collect().map { data ->
        def jsonContent = "[\n" + data.collect { it -> groovy.json.JsonOutput.toJson(it) }.join(",\n") + "\n]"
        def filePath = "${params.out_dir}/genome_info.json"
        new File(filePath).write(jsonContent)
        return filePath
    }

    emit:
    annotate_out.map {
        genome_name, gff, _tbl, _full_length_gff, _copies ->
        tuple(genome_name, gff)
    }

}