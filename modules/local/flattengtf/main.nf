process SUBREAD_FLATTENGTF {
    tag "${annotation}"
    label 'process_single'

    conda "bioconda::subread=2.0.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0'
        : 'biocontainers/subread:2.0.1--hed695b0_0'}"

    input:
    path annotation

    output:
    path "annotation.saf", emit: saf
    tuple val("${task.process}"), val('subread'), eval('flattenGTF -v 2>&1 | head -1 | sed -e "s/flattenGTF v//g"'), topic: versions, emit: versions_subread

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    flattenGTF ${args} -a ${annotation} -o annotation.saf
    """
}
