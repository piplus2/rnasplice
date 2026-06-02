process DEXSEQ_COUNT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f8/f894c18bfa323285c7a15ff518a8747ef5309d9c846afa81597688363c663dd3/data' :
        'community.wave.seqera.io/library/htseq:2.1.2--66d3f470b150ca54' }"

    input:
    tuple val(meta), path(bam), path (gff)
    val alignment_quality                   // val params.alignment_quality

    output:
    tuple val(meta), path("*.clean.count.txt"), emit: dexseq_clean_txt
    tuple val("${task.process}"), val('htseq'), eval('python -c "import HTSeq; print(HTSeq.__version__)"'), topic: versions, emit: versions_htseq

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def read_type = meta.single_end ? '' : '-p yes'

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = '-s yes'
    } else if (meta.strandedness == 'reverse') {
        strandedness = '-s reverse'
    } else if (meta.strandedness == 'unstranded') {
        strandedness = '-s no'
    }

    """
    dexseq_count.py ${gff} ${args} ${read_type} -f bam ${bam} -r pos ${prefix}.clean.count.txt -a ${alignment_quality} ${strandedness}
    """
}
