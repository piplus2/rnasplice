process DEXSEQ_ANNOTATION {
    tag "${gtf}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f8/f894c18bfa323285c7a15ff518a8747ef5309d9c846afa81597688363c663dd3/data' :
        'community.wave.seqera.io/library/htseq:2.1.2--66d3f470b150ca54' }"

    input:
    path gtf         // path gtf file
    val aggregation  // val params.aggregation

    output:
    path "*.gff"        , emit: gff
    tuple val("${task.process}"), val('htseq'), eval('python -c "import HTSeq; print(HTSeq.__version__)"'), topic: versions, emit: versions_htseq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "DEXSeq"
    def aggregation_arg = aggregation ? '' : '-r no'
    """
    dexseq_prepare_annotation.py ${gtf} ${prefix}.gff ${aggregation_arg} ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "DEXSeq"
    """
    echo $args

    touch ${prefix}.gff
    """
}
