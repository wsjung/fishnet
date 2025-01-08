nextflow.enable.dsl=2

include { RANDOM_PERMUTATION } from './modules.nf'
include { PREPROCESS_FOR_PASCAL } from './modules.nf'
include { RUN_PASCAL } from './modules.nf'
include { POSTPROCESS_PASCAL_OUTPUT } from './modules.nf'
include { GO_ANALYSIS } from './modules.nf'
include { MERGE_RESULTS } from './modules.nf'

workflow SUBWORKFLOW_MODULE_ENRICHMENT {

    take:
    gs
    module
    go

    main:
    RUN_PASCAL (
        gs | flatten,
        module | flatten,
        go | flatten
    )

    POSTPROCESS_PASCAL_OUTPUT (
        RUN_PASCAL.out.pascaloutput | flatten,
        RUN_PASCAL.out.genescorefile | flatten,
        RUN_PASCAL.out.gofile | flatten
    )

    GO_ANALYSIS (
        POSTPROCESS_PASCAL_OUTPUT.out.summaryslice | flatten,
        POSTPROCESS_PASCAL_OUTPUT.out.sigmodules | flatten,
        POSTPROCESS_PASCAL_OUTPUT.out.gofile | flatten
    )

    MERGE_RESULTS (
        GO_ANALYSIS.out.mastersummaryslice  | flatten,
        GO_ANALYSIS.out.gosummaries | flatten,
        GO_ANALYSIS.out.gofile | flatten
    )
}

workflow WORKFLOW_ORIGINAL_RUN {

    PREPROCESS_FOR_PASCAL (
        params.pvalFileName
    )

    SUBWORKFLOW_MODULE_ENRICHMENT (
        PREPROCESS_FOR_PASCAL.out.gs,
        PREPROCESS_FOR_PASCAL.out.module,
        PREPROCESS_FOR_PASCAL.out.go
    )

}

workflow WORKFLOW_RANDOM_RUN {

    RANDOM_PERMUTATION ()

    PREPROCESS_FOR_PASCAL (
        RANDOM_PERMUTATION.out.rp_scores | flatten
    )

    SUBWORKFLOW_MODULE_ENRICHMENT (
        PREPROCESS_FOR_PASCAL.out.gs,
        PREPROCESS_FOR_PASCAL.out.module,
        PREPROCESS_FOR_PASCAL.out.go
    )
}
