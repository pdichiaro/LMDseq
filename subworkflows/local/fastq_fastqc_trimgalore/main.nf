#!/usr/bin/env nextflow

//
// Cat fastq files and rname them according to the SHORT_NAME in sample sheet
//

include { CAT_FASTQ } from '../../../modules/local/cat/fastq/main'
include { FASTQC           } from '../../../modules/local/fastqc/main'
include { TRIMGALORE       } from '../../../modules/local/trimgalore/main'
include { FASTQC as FASTQC_TRIM } from '../../../modules/local/fastqc/main'

workflow FASTQ_FASTQC_TRIMGALORE {
    take:
    ch_input
    skip_fastqc
    skip_trimming

    main:

    ch_versions = Channel.empty()

    fastqc_html_raw = Channel.empty()
    fastqc_zip_raw  = Channel.empty()
    
    trim_unpaired   = Channel.empty()
    trim_reads      = Channel.empty()
    trim_html       = Channel.empty()
    trim_zip        = Channel.empty()
    trim_log        = Channel.empty()

    fastqc_html_trimmed = Channel.empty()
    fastqc_zip_trimmed  = Channel.empty()
        
    // concatenate fastqs
    CAT_FASTQ(ch_input)
    reads = CAT_FASTQ.out.reads
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    // read QC on raw
    if (!skip_fastqc) {
        FASTQC(reads)
        fastqc_html_raw = FASTQC.out.html
        fastqc_zip_raw  = FASTQC.out.zip
        ch_versions     = ch_versions.mix(FASTQC.out.versions.first())
    }

    // trimming  
    if (!skip_trimming) {
        TRIMGALORE(reads) 
        trim_unpaired = TRIMGALORE.out.unpaired
        trim_reads    = TRIMGALORE.out.reads
        trim_html     = TRIMGALORE.out.html
        trim_zip      = TRIMGALORE.out.zip
        trim_log      = TRIMGALORE.out.log
        ch_versions   = ch_versions.mix(TRIMGALORE.out.versions.first())
    }

    // read QC after trimming
    if (!skip_fastqc) {
        FASTQC_TRIM(trim_reads)
        fastqc_html_trimmed = FASTQC_TRIM.out.html
        fastqc_zip_trimmed  = FASTQC_TRIM.out.zip
        ch_versions         = ch_versions.mix(FASTQC_TRIM.out.versions.first())
    }

    emit:
    reads
    trim_reads

    raw_fastqc_html = fastqc_html_raw
    raw_fastqc_zip  = fastqc_zip_raw

    trim_fastqc_html = fastqc_html_trimmed
    trim_fastqc_zip  = fastqc_zip_trimmed

    trim_unpaired
    trim_html
    trim_zip
    trim_log

    versions = ch_versions
    multiqc_files = fastqc_html_raw.mix(fastqc_zip_raw)
                                  .mix(fastqc_html_trimmed)
                                  .mix(fastqc_zip_trimmed)
                                  .mix(trim_html)
                                  .mix(trim_zip)
                                  .mix(trim_log)
}