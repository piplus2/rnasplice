process DEXSEQ_ANNOTATION {
    tag "$gtf"
    label 'process_medium'

    conda "bioconda::htseq=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/htseq:2.0.2--py310ha14a713_0' :
    'biocontainers/htseq:2.0.2--py310ha14a713_0' }"

    input:
    path gtf         // path gtf file
    val aggregation  // val params.aggregation

    output:
    path "*.gff"        , emit: gff
    tuple val("${task.process}"), val('htseq'), eval('python -c "import pkg_resources; print(pkg_resources.get_distribution(\"htseq\").version)"'), topic: versions, emit: versions_htseq

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "DEXSeq"

    def aggregation_arg = aggregation ? '' : '-r no'

    """
    dexseq_prepare_annotation.py ${gtf} ${prefix}.gff ${aggregation_arg} ${args}
    """
}
