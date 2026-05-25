process LEAFCUTTER_CLUSTER {
    label 'process_medium'

    conda "conda-forge::python=3.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.7' :
        'quay.io/biocontainers/python:3.7' }"

    input:
    path junc_files
    path gtf

    output:
    path("*perind_numers.counts.gz") , emit: counts
    path("*perind.counts.gz")        , emit: perind_counts
    tuple val("${task.process}"), val('leafcutter'), eval('leafcutter_cluster_regtools.py --version 2>&1 | head -n 1 | sed "s/leafcutter_cluster_regtools.py v//g"'), topic: versions, emit: versions_leafcutter


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ls *junc > test_juncfiles.txt
    leafcutter_cluster_regtools.py \\
        $args \\
        -j test_juncfiles.txt \\
        -o lc
    """
}
