#!/usr/bin/env nextflow

include { RNASEQ } from './workflows/rnaseq.nf'
include { BRAKER3 } from './workflows/braker3.nf'
include { INTERPRO } from './workflows/interpro.nf'
include { GQUALS } from './workflows/gquals.nf'
include { GRAPHMAP } from './workflows/graphmap.nf'
include { GRAPHCALL_GAM } from './workflows/graphcall.nf'
include { GRAPHCALL_BCFMERGE } from './workflows/graphcall.nf'
include { NODE_CALL } from './workflows/gwas.nf'
include { ASSEMPREP } from './workflows/assemprep.nf'

def helpMessage() {

    log.info"""

    ccdm-bnpan - NEXTFLOW WORKFLOW version ${workflow.manifest.version}
    =================================

    Parameters:

        Required:

        --workflow              One of either:
                                    'rnaseq' 'braker3' 'interproscan' 'gquals'
                                    'graphmap' 'graphcallg' 'bcfmerge' 'gwas'
                                    'assemprep'

        Optional:

        --help                  Show this help message and exit.
        --version               Show the version of the pipeline and exit.
    
    """.stripIndent()
}

def error_check() {

    // Check if a valid workflow is specified.
    def valid_workflows = 
        [
            'rnaseq', 'braker3', 'interproscan',
            'gquals', 'graphmap', 'graphcallg',
            'bcfmerge', 'gwas', 'assemprep'
         ]

    if (!params.workflow) {

        if (params.help) {

            helpMessage()
            System.exit(0)

        } else {

            error "Workflow ${params.workflow} is not recognised."
        }

    } else if (!valid_workflows.contains(params.workflow)) {

        error "Workflow ${params.workflow} is not recognised."

    }

}

workflow {

    error_check()

    if (params.workflow == 'rnaseq') {

        RNASEQ()

    } else if (params.workflow == 'braker3') {

        BRAKER3()

    } else if (params.workflow == 'interproscan') {

        INTERPRO()

    } else if (params.workflow == 'gquals') {

        GQUALS()

    } else if (params.workflow == 'graphmap') {

        GRAPHMAP()
    
    } else if (params.workflow == 'graphcallg') {

        GRAPHCALL_GAM()
    
    } else if (params.workflow == 'bcfmerge') {

        GRAPHCALL_BCFMERGE()

    } else if (params.workflow == 'gwas') {

        NODE_CALL()

    } else if (params.workflow == 'assemprep') {

        ASSEMPREP()

    } else {

        error "Workflow ${params.workflow} is not recognised. Valid options are: rnaseq, braker3, interproscan, gquals, graphmap, graphcallg, bcfmerge, gwas, assemprep."

    }

}