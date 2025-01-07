#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { PREPROCESS_FOR_PASCAL } from './modules.nf'

workflow {
    PREPROCESS_FOR_PASCAL(params.pvalFileName)
}
