process STAGER {
    label 'process_medium'

    conda "bioconda::bioconductor-stager=1.12.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-stager:1.12.0--r40_0' :
        'biocontainers/bioconductor-stager:1.12.0--r40_0' }"

    input:
    tuple val(contrast), path(feature_tsv), path(gene_tsv)
    val analysis_type                                         // val: "dexseq" or "drimseq"

    output:
    path "stageRTx.*.rds"              , emit: stager_rds
    path "getAdjustedPValues.*.rds"    , emit: stager_padj_rds
    path "getAdjustedPValues.*.tsv"    , emit: stager_padj_tsv
    tuple val("${task.process}"), val('stageR'), eval('Rscript -e "library(stageR); cat(as.character(packageVersion(\'stageR\')))"'), topic: versions, emit: versions_stager

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    run_stager.R ${contrast} ${feature_tsv} ${gene_tsv} ${analysis_type} ${args}
    """
}
