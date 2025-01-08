#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { WORKFLOW_ORIGINAL_RUN } from './subworkflows.nf'
include { WORKFLOW_RANDOM_RUN } from './subworkflows.nf'

workflow {
    if (params.random_permutation) {
        WORKFLOW_RANDOM_RUN()
    } else {
        WORKFLOW_ORIGINAL_RUN()
    }
}
