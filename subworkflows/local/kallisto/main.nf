#!/usr/bin/env nextflow

//
// Pseudo-alignment and quantification by Kallisto
//

include { KALLISTO_QUANT      } from '../../../modules/local/kallisto/quant'
include { BAM_COVERAGE        } from '../../../modules/local/bam_coverage/main'
include { SAMTOOLS            } from '../../../modules/local/samtools/main'
include { TXIMPORT            } from '../../../modules/local/tximport/main'
include { DESEQ2_NORM_QC      } from '../../../modules/local/deseq2_norm_qc/main'


workflow KALLISTO {
    take:
    reads                        // channel: [ val(meta), [ reads ] ]
    index                        // channel: /path/to//index/
    gtf                          // channel: /path/to/genome.gtf
    chromosomes                  // channel: /path/to/Chr.txt
    kallisto_quant_fraglen       // val: Estimated fragment length required by Kallisto in single-end mode
    kallisto_quant_fraglen_sd    // val: Estimated standard error for fragment length required by Kallisto in single-end mode
    bin_size                     // val: 
    reference                    // channel: /path/to/hg38_ref.txt
    sample_metadata              // channel: /path/to/Sample.txt

    main:

    ch_versions = Channel.empty()
    fastqc_html = Channel.empty()
    ch_pseudo_results = Channel.empty()
    ch_bigwig = Channel.empty()
    ch_matrix = Channel.empty()
    ch_qc_files = Channel.empty()

        
    // Map reads with kallisto
    KALLISTO_QUANT (
        reads,
        index,
        gtf,
        chromosomes,
        kallisto_quant_fraglen,
        kallisto_quant_fraglen_sd
    )
    ch_pseudo_results = KALLISTO_QUANT.out.results
    ch_pseudo_bam = KALLISTO_QUANT.out.bam
    ch_pseudo_bai = KALLISTO_QUANT.out.bai
    ch_pseudo_multiqc = KALLISTO_QUANT.out.log
    ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions.first())


    // Define a scale factor on the mapped reads
    SAMTOOLS (
        ch_pseudo_bam
    )
    ch_scale_f = SAMTOOLS.out.scale_f.map { meta, stdout_val -> tuple(meta, stdout_val.trim())}    
    ch_bam_with_bai = ch_pseudo_bam.join(ch_pseudo_bai)
    ch_bam_cov = ch_bam_with_bai.join(ch_scale_f)
    ch_versions = ch_versions.mix(SAMTOOLS.out.versions)


    // Create bigwig scaled to the defined scale factor
    BAM_COVERAGE (
        ch_bam_cov,
        bin_size
    )
    ch_bigwig = BAM_COVERAGE.out.bigwig
    ch_versions = ch_versions.mix(BAM_COVERAGE.out.versions)

    ch_all_quant = ch_pseudo_results.map { it[1] }.collect()  // Collect all results


    // Create a reference database 
    TXIMPORT (
        ch_all_quant,
        gtf,
        reference
    )
    ch_matrix = TXIMPORT.out.matrix
    ch_versions = ch_versions.mix(TXIMPORT.out.versions)

    // DESeq2 normalization and quality control
    DESEQ2_NORM_QC (
        ch_matrix,
        sample_metadata
    )
    ch_norm_files = DESEQ2_NORM_QC.out.norm_files
    ch_count_files = DESEQ2_NORM_QC.out.count_files
    ch_qc_files = DESEQ2_NORM_QC.out.qc_files
    ch_normalized_counts = DESEQ2_NORM_QC.out.normalized_counts
    ch_rlog_counts = DESEQ2_NORM_QC.out.rlog_counts
    ch_dds_rdata = DESEQ2_NORM_QC.out.dds_rdata
    ch_versions = ch_versions.mix(DESEQ2_NORM_QC.out.versions)


    emit:

    results                       = ch_pseudo_results    // channel: [ val(meta), results_dir ]

    bigwig                        = ch_bigwig            // channel: [ val(meta), scaled_bigwig ]
    matrix                        = ch_matrix            // channel: [ val(meta), gene_expression_matrix ]
    norm_files                    = ch_norm_files        // channel: [ path(6_Norm_folder/**) ]
    count_files                   = ch_count_files       // channel: [ path(7_Counts_folder/**) ]
    qc_files                      = ch_qc_files          // channel: [ path(8_Quality_folder/**) ]
    normalized_counts             = ch_normalized_counts // channel: [ path(deseq2_normalized_counts.txt) ]
    rlog_counts                   = ch_rlog_counts       // channel: [ path(deseq2_rlog_counts.txt) ]
    dds_rdata                     = ch_dds_rdata         // channel: [ path(deseq2.dds.RData) ]

    versions                      = ch_versions          // channel: [ versions.yml ]
    multiqc_files                 = ch_pseudo_multiqc    // channel: [ val(meta), files_for_multiqc ]

}

