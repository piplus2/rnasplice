//
// SUPPA Subworkflow
//

include { GENERATE_EVENTS as GENERATE_EVENTS_IOE } from '../../../modules/local/suppa_generateevents'
include { GENERATE_EVENTS as GENERATE_EVENTS_IOI } from '../../../modules/local/suppa_generateevents'

include { PSIPEREVENT   } from '../../../modules/local/suppa_psiperevent'
include { PSIPERISOFORM } from '../../../modules/local/suppa_psiperisoform'

include { SPLIT_FILES as SPLIT_FILES_TPM } from '../../../modules/local/suppa_split_files'
include { SPLIT_FILES as SPLIT_FILES_IOE } from '../../../modules/local/suppa_split_files'
include { SPLIT_FILES as SPLIT_FILES_IOI } from '../../../modules/local/suppa_split_files'

include { DIFFSPLICE as DIFFSPLICE_IOE } from '../../../modules/local/suppa_diffsplice'
include { DIFFSPLICE as DIFFSPLICE_IOI } from '../../../modules/local/suppa_diffsplice'

include { CLUSTERGROUPS as CLUSTERGROUPS_IOE } from '../../../modules/local/suppa_clustergroups'
include { CLUSTERGROUPS as CLUSTERGROUPS_IOI } from '../../../modules/local/suppa_clustergroups'

include { CLUSTEREVENTS as CLUSTEREVENTS_IOE } from '../../../modules/local/suppa_clusterevents'
include { CLUSTEREVENTS as CLUSTEREVENTS_IOI } from '../../../modules/local/suppa_clusterevents'

workflow SUPPA {

    take:

    ch_gtf
    ch_tpm
    ch_samplesheet
    ch_contrastsheet
    suppa_per_local_event       // params.suppa_per_local_event
    generateevents_boundary     // params.generateevents_boundary
    generateevents_threshold    // params.generateevents_threshold
    generateevents_exon_length  // params.generateevents_exon_length
    generateevents_event_type   // params.generateevents_event_type
    generateevents_pool_genes   // params.generateevents_pool_genes
    psiperevent_total_filter    // params.psiperevent_total_filter
    diffsplice_local_event      // params.diffsplice_local_event
    diffsplice_method           // params.diffsplice_method
    diffsplice_area             // params.diffsplice_area
    diffsplice_lower_bound      // params.diffsplice_lower_bound
    diffsplice_alpha            // params.diffsplice_alpha
    diffsplice_tpm_threshold    // params.diffsplice_tpm_threshold
    diffsplice_nan_threshold    // params.diffsplice_nan_threshold
    diffsplice_gene_correction  // params.diffsplice_gene_correction
    diffsplice_paired           // params.diffsplice_paired
    diffsplice_median           // params.diffsplice_median
    clusterevents_local_event   // params.clusterevents_local_event
    clusterevents_dpsithreshold // params.clusterevents_dpsithreshold
    clusterevents_eps           // params.clusterevents_eps
    clusterevents_metric        // params.clusterevents_metric
    clusterevents_min_pts       // params.clusterevents_min_pts
    clusterevents_method        // params.clusterevents_method
    clusterevents_sigthreshold  // params.clusterevents_sigthreshold
    clusterevents_separation    // params.clusterevents_separation
    suppa_per_isoform           // params.suppa_per_isoform

    main:

    // Split the tpm file (contains all samples) into individual files based on condition

    def output_type = ".tpm"
    def calc_ranges = false
    def prefix = ''
    def file_type = ''

    SPLIT_FILES_TPM (
        ch_tpm,
        ch_samplesheet,
        output_type,
        calc_ranges,
        prefix
    )

    ch_split_suppa_tpms = SPLIT_FILES_TPM.out.tpms

    // If per AS local analysis:

    ch_ioe_events             = channel.empty()
    ch_suppa_local_psi        = channel.empty()
    ch_split_suppa_local_psi  = channel.empty()

    ch_dpsi_local             = channel.empty()
    ch_psivec_local           = channel.empty()

    ch_groups_ioe             = channel.empty()
    ch_cluster_vec_local      = channel.empty()
    ch_cluster_log_local      = channel.empty()

    if (suppa_per_local_event) {

        file_type = 'ioe'

        // Generate AS events on the GTF - Local events

        GENERATE_EVENTS_IOE (
            ch_gtf,
            file_type,
            generateevents_boundary,
            generateevents_threshold,
            generateevents_exon_length,
            generateevents_event_type,
            generateevents_pool_genes
        )

        ch_ioe_events  = GENERATE_EVENTS_IOE.out.events

        // Calculate the psi values of Local events (using events file and TPM)

        PSIPEREVENT (
            ch_ioe_events,
            ch_tpm,
            psiperevent_total_filter
        )

        ch_suppa_local_psi = PSIPEREVENT.out.psi

        // Split the PSI files between the conditions

        output_type = ".psi"
        calc_ranges = true
        prefix = "local"

        SPLIT_FILES_IOE (
            ch_suppa_local_psi,
            ch_samplesheet,
            output_type,
            calc_ranges,
            prefix
        )

        ch_split_suppa_local_psi = SPLIT_FILES_IOE.out.psis

        // Calculate differential analysis between conditions

        if (diffsplice_local_event) {

            // Create contrasts channel

            ch_suppa_local_contrasts_raw = ch_contrastsheet

            // Add TPM files to contrasts channel

            ch_suppa_tpm_conditions = SPLIT_FILES_TPM.out.tpms
                .flatten()
                .map { [it.baseName, it ] }

            ch_suppa_local_contrasts_tpm1 = ch_suppa_local_contrasts_raw
                .map { it -> [it['treatment'], it] }
                .combine ( ch_suppa_tpm_conditions, by: 0 )
                .map { it -> it[1] + ['tpm1': it[2]] }

            ch_suppa_local_contrasts_tpm2 = ch_suppa_local_contrasts_tpm1
                .map { it -> [it['control'], it] }
                .combine ( ch_suppa_tpm_conditions, by: 0 )
                .map { it -> it[1] + ['tpm2': it[2]] }

            // Add PSI files to contrasts channel

            ch_suppa_psi_conditions = SPLIT_FILES_IOE.out.psis
                .flatten()
                .map { [ it.baseName.toString().replaceAll("local_", ""), it ] }

            ch_suppa_local_contrasts_psi1 = ch_suppa_local_contrasts_tpm2
                .map { it -> [it['treatment'], it] }
                .combine ( ch_suppa_psi_conditions, by: 0 )
                .map { it -> it[1] + ['psi1': it[2]] }

            ch_suppa_local_contrasts_psi2 = ch_suppa_local_contrasts_psi1
                .map { it -> [it['control'], it] }
                .combine ( ch_suppa_psi_conditions, by: 0 )
                .map { it -> it[1] + ['psi2': it[2]] }
                .map { it -> [it['control'], it] }
                .combine ( ch_suppa_psi_conditions, by: 0 )
                .map { it -> it[1] + ['psi2': it[2]] }

            // Create input channels to diffsplice process

            ch_split_suppa_tpms_psi = ch_suppa_local_contrasts_psi2.map { it -> [ it.treatment, it.control, it.tpm1, it.tpm2, it.psi1, it.psi2 ] }

            DIFFSPLICE_IOE(
                ch_ioe_events,
                ch_split_suppa_tpms_psi,
                prefix,
                diffsplice_method,
                diffsplice_area,
                diffsplice_lower_bound,
                diffsplice_alpha,
                diffsplice_tpm_threshold,
                diffsplice_nan_threshold,
                diffsplice_gene_correction,
                diffsplice_paired,
                diffsplice_median
            )

            ch_dpsi_local     = DIFFSPLICE_IOE.out.dpsi
            ch_psivec_local   = DIFFSPLICE_IOE.out.psivec

            if (clusterevents_local_event) {

                // Get ranges for cluster analysis

                CLUSTERGROUPS_IOE ( ch_psivec_local )

                ch_groups_ioe = CLUSTERGROUPS_IOE.out[0]

                // Join channels to ensure consistent order

                ch_clusterevents_ioe = ch_dpsi_local.join(ch_psivec_local, by: [0, 1]).join(ch_groups_ioe, by: [0, 1])

                // Run Clustering

                CLUSTEREVENTS_IOE(
                    ch_clusterevents_ioe,
                    prefix,
                    clusterevents_dpsithreshold,
                    clusterevents_eps,
                    clusterevents_metric,
                    clusterevents_min_pts,
                    clusterevents_method,
                    clusterevents_sigthreshold,
                    clusterevents_separation
                )

                ch_cluster_vec_local   = CLUSTEREVENTS_IOE.out.clustvec
                ch_cluster_log_local   = CLUSTEREVENTS_IOE.out.cluster_log
            }
        }
    }

    // If per isoform analysis:

    ch_ioi_events              = channel.empty()
    ch_suppa_isoform_psi       = channel.empty()
    ch_split_suppa_isoform_psi = channel.empty()

    ch_dpsi_isoform            = channel.empty()
    ch_psivec_isoform          = channel.empty()

    ch_groups_ioi              = channel.empty()
    ch_cluster_vec_isoform     = channel.empty()
    ch_cluster_log_isoform     = channel.empty()

    if (suppa_per_isoform) {

        file_type = 'ioi'

        // Generate events - transcript level

        GENERATE_EVENTS_IOI (
            ch_gtf,
            file_type,
            generateevents_boundary,
            generateevents_threshold,
            generateevents_exon_length,
            generateevents_event_type,
            generateevents_pool_genes
        )

        ch_ioi_events = GENERATE_EVENTS_IOI.out.events

        // Get the psi values per isoform

        PSIPERISOFORM (
            ch_gtf,
            ch_tpm
        )

        ch_suppa_isoform_psi = PSIPERISOFORM.out.psi

        // Split the PSI files between the conditions

        output_type = ".psi"
        calc_ranges = true
        prefix = "transcript"

        SPLIT_FILES_IOI (
            ch_suppa_isoform_psi,
            ch_samplesheet,
            output_type,
            calc_ranges,
            prefix
        )

        ch_split_suppa_isoform_psi = SPLIT_FILES_IOI.out.psis

        // Calculate differential analysis between conditions - Transcript level

        if (params.diffsplice_isoform) {

            // Create contrasts channel

            ch_suppa_isoform_contrasts_raw = ch_contrastsheet

            // Add TPM files to contrasts channel

            ch_suppa_tpm_conditions = SPLIT_FILES_TPM.out.tpms
                .flatten()
                .map { [it.baseName, it ] }

            ch_suppa_isoform_contrasts_tpm1 = ch_suppa_isoform_contrasts_raw
                .map { it -> [it['treatment'], it] }
                .combine ( ch_suppa_tpm_conditions, by: 0)
                .map { it -> it[1] + ['tpm1': it[2]] }

            ch_suppa_isoform_contrasts_tpm2 = ch_suppa_isoform_contrasts_tpm1
                .map { it -> [it['control'], it] }
                .combine ( ch_suppa_tpm_conditions, by: 0)
                .map { it -> it[1] + ['tpm2': it[2]] }

            // Add PSI files to contrasts channel

            ch_suppa_psi_conditions = SPLIT_FILES_IOI.out.psis
                .flatten()
                .map { [ it.baseName.toString().replaceAll("transcript_", ""), it ] }

            ch_suppa_isoform_contrasts_psi1 = ch_suppa_isoform_contrasts_tpm2
                .map { it -> [it['treatment'], it] }
                .combine ( ch_suppa_psi_conditions, by: 0 )
                .map { it -> it[1] + ['psi1': it[2]] }

            ch_suppa_isoform_contrasts_psi2 = ch_suppa_isoform_contrasts_psi1
                .map { it -> [it['control'], it] }
                .combine ( ch_suppa_psi_conditions, by: 0 )
                .map { it -> it[1] + ['psi2': it[2]] }

            // Create input channels to diffsplice process

            ch_split_suppa_tpms = ch_suppa_isoform_contrasts_psi2.map { it -> [ it.treatment, it.control, it.tpm1, it.tpm2, it.psi1, it.psi2 ] }

            DIFFSPLICE_IOI(
                ch_ioi_events,
                ch_split_suppa_tpms_psi,
                prefix,
                diffsplice_method,
                diffsplice_area,
                diffsplice_lower_bound,
                diffsplice_alpha,
                diffsplice_tpm_threshold,
                diffsplice_nan_threshold,
                diffsplice_gene_correction,
                diffsplice_paired,
                diffsplice_median
            )

            ch_dpsi_isoform   = DIFFSPLICE_IOI.out.dpsi
            ch_psivec_isoform = DIFFSPLICE_IOI.out.psivec

            if (params.clusterevents_isoform) {

                // Get ranges for cluster analysis

                CLUSTERGROUPS_IOI ( ch_psivec_isoform )

                ch_groups_ioi = CLUSTERGROUPS_IOI.out[0]

                // Join channels to ensure consistent order

                ch_clusterevents_ioi = ch_dpsi_isoform.join(ch_psivec_isoform, by: [0, 1]).join(ch_groups_ioi, by: [0, 1])

                // Run Clustering

                CLUSTEREVENTS_IOI(
                    ch_clusterevents_ioi,
                    prefix,
                    clusterevents_dpsithreshold,
                    clusterevents_eps,
                    clusterevents_metric,
                    clusterevents_min_pts,
                    clusterevents_method,
                    clusterevents_sigthreshold,
                    clusterevents_separation
                )

                ch_cluster_vec_isoform = CLUSTEREVENTS_IOI.out.clustvec
                ch_cluster_log_isoform = CLUSTEREVENTS_IOI.out.cluster_log

            }
        }
    }

    // Define output

    emit:

    ioe_events              = ch_ioe_events                      //    path: events ioe
    ioi_events              = ch_ioi_events                      //    path: events ioi

    suppa_local_psi         = ch_suppa_local_psi                 //    path: suppa_local.psi
    suppa_isoform_psi       = ch_suppa_isoform_psi               //    path: suppa_isoform.psi

    split_suppa_tpms        = ch_split_suppa_tpms                //    path: suppa_cond1.tpm, suppa_cond2.tpm
    split_suppa_local_psi   = ch_split_suppa_local_psi           //    path: suppa_local_cond1.psi, suppa_local_cond2.psi
    split_suppa_isoform_psi = ch_split_suppa_isoform_psi         //    path: suppa_isoform_cond1.psi, suppa_isoform_cond2.psi

    dpsi_local              = ch_dpsi_local                      //    path: local.dpsi
    psivec_local            = ch_psivec_local                    //    path: local.psivec
    dpsi_isoform            = ch_dpsi_isoform                    //    path: isoform.dpsi
    psivec_isoform          = ch_psivec_isoform                  //    path: isoform.psivec

    cluster_vec_local       = ch_cluster_vec_local               //    path: local.clustvec
    cluster_log_local       = ch_cluster_log_local               //    path: local.log
    cluster_vec_isoform     = ch_cluster_vec_isoform             //    path: isoform.clustvec
    cluster_log_isoform     = ch_cluster_log_isoform             //    path: isoform.log
}
