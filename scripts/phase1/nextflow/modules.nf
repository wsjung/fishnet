nextflow.enable.dsl=2

process PREPROCESS_FOR_PASCAL {

    //conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9d836da785124bb367cbe6fbfc00dddd2107a4da:b033d6a4ea3a42a6f5121a82b262800f1219b382-0' :
        'quay.io/biocontainers/mulled-v2-9d836da785124bb367cbe6fbfc00dddd2107a4da:b033d6a4ea3a42a6f5121a82b262800f1219b382-0' }"

    //label "process_low"
    publishDir "./results/", pattern: "pascalInput/*", mode: 'copy'

    input:
    path pvalFile

    output:
    path("pascalInput/GS_*")
    path("pascalInput/Module_*")
    path("pascalInput/GO_*")

    script:
    """
    python3 ${projectDir}/bin/preProcessForPascal.py \
        ${pvalFile} \
        ${params.moduleFileDir} \
        "pascalInput/" \
        ${params.pipeline} \
        ${params.trait} \
        ${params.geneColName} \
        ${params.pvalColName}
    """
}
