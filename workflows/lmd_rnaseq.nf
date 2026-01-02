#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
// MultiQC removed - not needed for LMDseq workflow


//
// SUBWORKFLOW: Consisting of local/modules
//
include { FASTQ_FASTQC_TRIMGALORE } from '../subworkflows/local/fastq_fastqc_trimgalore/main' 
include { KALLISTO } from '../subworkflows/local/kallisto/main' 


//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { paramsSummaryMap                 } from 'plugin/nf-schema'
include { samplesheetToList                } from 'plugin/nf-schema'
include { softwareVersionsToYAML           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { checkSamplesAfterGrouping        } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { multiqcTsvFromList               } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { biotypeInGtf                     } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Header files for MultiQC
//ch_pca_header_multiqc           = file("$projectDir/workflows/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
//sample_status_header_multiqc    = file("$projectDir/workflows/assets/multiqc/sample_status_header.txt", checkIfExists: true)
//ch_clustering_header_multiqc    = file("$projectDir/workflows/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
//ch_biotypes_header_multiqc      = file("$projectDir/workflows/assets/multiqc/biotypes_header.txt", checkIfExists: true)
//ch_dummy_file                   = ch_pca_header_multiqc

workflow RNASEQ {
    
    take:
    //ch_samplesheet        // channel: path(sample_sheet.txt)
    ch_versions             // channel: [ path(versions.yml) ]
    ch_fasta                // channel: path(genome.fasta)
    ch_gtf                  // channel: path(gtf)
    ch_transcript_fasta     // channel: path(transcript.fasta)
    ch_index                // channel: [ meta, path(kallisto/index/) ]
    ch_chrom_sizes          // channel: path(genome.chrom.sizes)

    main:
    ch_versions = Channel.empty()

    // -----------------------
    // Run FASTQ preprocessing
    // -----------------------

    // Create channel from input file provided through ch_samplesheet
    ch_samplesheet = Channel.fromPath(params.input)
    ch_input = ch_samplesheet
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            def fastq_dir = file(row.FASTQ_FOLDER)
            def is_single_end = params.seq_mode == 'SE'

            def meta = [
                //id: row.SID,
                id: "${row.SHORT_NAME}_${row.REPLICATE}",
                short_name: row.SHORT_NAME,
                replicate: row.REPLICATE,
                single_end: is_single_end,
                strandedness: row.strandedness
            ]

            // Create file objects and handle file discovery robustly
            def files = []
            
            // Try to find files in the directory
            if (fastq_dir.exists() && fastq_dir.isDirectory()) {
                def all_files = fastq_dir.listFiles()?.findAll { it.name.endsWith('.fastq.gz') } ?: []
                def r1_files = all_files.findAll { it.name.contains('_R1_') }
                def r2_files = all_files.findAll { it.name.contains('_R2_') }
                
                if (is_single_end && r1_files) {
                    files = r1_files
                } else if (!is_single_end && r1_files && r2_files) {
                    files = r1_files + r2_files
                } else if (!is_single_end && r1_files) {
                    // PE mode but only R1 found, use R1 only
                    files = r1_files
                }
            }
            
            // Skip this sample if no valid files found
            if (!files || files.isEmpty()) {
                return null
            }
            
            [ meta, files ]
        }
        .filter { it != null }
    
    // No need to fork channel anymore since MultiQC has been removed


    // -----------------------
    // Run Subworkflow: cat FASTQs, rename files, read QC and trim adapters
    // -----------------------

    if (params.trimmer == 'trimgalore') {
        FASTQ_FASTQC_TRIMGALORE (
            ch_input,
            params.skip_fastqc,
            params.skip_trimming
        )
    }
    ch_versions = ch_versions.mix(FASTQ_FASTQC_TRIMGALORE.out.versions)
    ch_trimmed_reads = FASTQ_FASTQC_TRIMGALORE.out.reads
    
    //ch_trimmed_reads.collect { println it } // This will give you the current contents of the channel


    // -----------------------
    // Run Subworkflow: Alignment with KALLISTO, generate pseudo-bam, raw count matrix, and DESeq2 normalization/QC
    // -----------------------
    if (params.pseudo_aligner == 'kallisto') {
        KALLISTO(
            ch_trimmed_reads,
            ch_index,
            ch_gtf,
            ch_chrom_sizes,
            params.kallisto_quant_fraglen,
            params.kallisto_quant_fraglen_sd,
            params.bin_size,
            params.reference ? Channel.value(file(params.reference)) : Channel.value([]),
            ch_samplesheet
        )
    }
    ch_versions = ch_versions.mix(KALLISTO.out.versions)
    ch_bigwig = KALLISTO.out.bigwig
    ch_raw_matrix = KALLISTO.out.matrix
    ch_norm_files = KALLISTO.out.norm_files
    ch_count_files = KALLISTO.out.count_files
    ch_qc_files = KALLISTO.out.qc_files
    ch_normalized_counts = KALLISTO.out.normalized_counts
    ch_rlog_counts = KALLISTO.out.rlog_counts
    ch_dds_rdata = KALLISTO.out.dds_rdata


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'LMD_rnaseq_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MultiQC has been removed from this workflow
    // All QC outputs are available in individual result folders
    //


    emit:
    raw_matrix          = ch_raw_matrix          // channel: [ path(EX_reads_RAW.txt) ]
    norm_files          = ch_norm_files          // channel: [ path(6_Norm_folder/**) ]
    count_files         = ch_count_files         // channel: [ path(7_Counts_folder/**) ]
    qc_files            = ch_qc_files            // channel: [ path(8_Quality_folder/**) ]
    versions            = ch_versions            // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
