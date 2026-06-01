process STAGER {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/41/41a73be6620eaa0bb67220f2b3e71168e4724c09556c76b65c2ad88f383d7bd4/data' :
        'community.wave.seqera.io/library/bioconductor-stager:1.32.0--854b4a08a81d122d' }"

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
