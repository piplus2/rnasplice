process GENERATE_EVENTS {
    tag "$gtf"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d887a6a05dec2a1f64fdff0eac40581f9a1ec30301b2c267bde7f564b0f14270/data' :
        'community.wave.seqera.io/library/suppa:2.4--2612fcca3884f6bc' }"

    input:
    path gtf
    val file_type
    val generateevents_boundary       // val params.generateevents_boundary
    val generateevents_threshold      // val params.generateevents_threshold
    val generateevents_exon_length    // val params.generateevents_exon_length
    val generateevents_event_type     // val params.generateevents_event_type
    val generateevents_pool_genes     // val params.generateevents_pool_genes

    output:
    path "events.*"     , emit: events
    path "events_*.*"   , emit: eventstype, optional : true
    tuple val("${task.process}"), val('suppa'), eval("suppa.py -v | sed '1!d;s/.* //'"), topic: versions, emit: versions_suppa

    when:
    task.ext.when == null || task.ext.when

    script:

    // If pool_genes is set to true then include the -p parameter
    def poolgenes = generateevents_pool_genes ? "-p" : ''

    // Calculate local events and combine all the ioe events files
    if (file_type == 'ioe') {
        """
        suppa.py \\
            generateEvents \\
            -i $gtf \\
            -f $file_type \\
            -o events \\
            -e $generateevents_event_type \\
            -b $generateevents_boundary \\
            -t $generateevents_threshold \\
            -l $generateevents_exon_length \\
            $poolgenes

        awk 'FNR==1 && NR!=1 { while (/^seqname/) getline; }  1 {print}' *.ioe > events.ioe
        """
    }
    // Calculate transcript events
    else if (file_type == 'ioi') {
        """
        suppa.py \\
            generateEvents \\
            -i $gtf \\
            -f $file_type \\
            -o events \\
            $poolgenes
        """
    }
}
