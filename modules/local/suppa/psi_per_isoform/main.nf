process PSIPERISOFORM {
    tag "$tpm"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d887a6a05dec2a1f64fdff0eac40581f9a1ec30301b2c267bde7f564b0f14270/data' :
        'community.wave.seqera.io/library/suppa:2.4--2612fcca3884f6bc' }"

    input:
    path gtf
    path tpm

    output:
    path "suppa_isoform.psi"  , emit: psi
    tuple val("${task.process}"), val('suppa'), eval("suppa.py -v | sed '1!d;s/.* //'"), topic: versions, emit: versions_suppa

    when:
    task.ext.when == null || task.ext.when

    script: // Calculate the psi values isoform level
    """
    suppa.py \\
        psiPerIsoform \\
        -g $gtf \\
        -e $tpm \\
        -o suppa
    """
}
