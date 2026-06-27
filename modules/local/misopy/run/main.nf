process MISOPY_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/misopy:0.5.4--py27h516909a_2':
        'quay.io/biocontainers/misopy:0.5.4--py27h516909a_2' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path miso_index

    output:
    tuple val(meta), path("miso_data/*")    , emit: miso
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //g"'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('misopy'), eval('python -c "import pkg_resources; print(pkg_resources.get_distribution(\'misopy\').version)"'), topic: versions, emit: versions_misopy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    miso --run \\
        ${miso_index} \\
        ${bam} \\
        --output-dir miso_data/${prefix} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir -p miso_data/${prefix}
    touch miso_data/${prefix}/miso_results
    """
}
