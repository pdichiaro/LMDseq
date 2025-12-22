// Simplified PREPARE_GENOME workflow with only KALLISTO as aligner

include { GUNZIP as GUNZIP_GTF }              from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_FASTA }            from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../../modules/nf-core/gunzip'
include { GFFREAD }                           from '../../../modules/nf-core/gffread'
include { GTF_FILTER }                        from '../../../modules/local/gtf_filter'
include { KALLISTO_INDEX }                    from '../../../modules/local/kallisto/index'
include { UNTAR as UNTAR_KALLISTO_INDEX }     from '../../../modules/nf-core/untar'
include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE } from '../../../modules/local/preprocess_transcripts_fasta_gencode'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA } from '../../../modules/nf-core/rsem/preparereference'
include { CUSTOM_GETCHROMSIZES }              from '../../../modules/local/custom_getchromsizes'

workflow PREPARE_GENOME {

    take:
    fasta                   // file: /path/to/genome.fasta (optional!)
    gtf                     // file: /path/to/genome.gtf
    gff                     // file: /path/to/genome.gff
    transcript_fasta        // file: /path/to/transcript.fasta
    index                   // directory: /path/to/kallisto/index/
    gencode                 // boolean: whether the genome is from GENCODE
    pseudo_aligner          // string: Specifies the pseudo-alignment algorithm to use - available option is kallisto
    skip_gtf_filter         // boolean: Skip filtering of GTF for valid scaffolds and/ or transcript IDs

    main:
    ch_versions = Channel.empty()

    // Handle GTF or GFF to GTF
    ch_gtf = Channel.empty()
    if (gtf) {
        ch_gtf = gtf.endsWith('.gz') 
            ? GUNZIP_GTF([ [:], file(gtf, checkIfExists: true) ]).gunzip.map { it[1] }
            : Channel.value(file(gtf, checkIfExists: true))
    } else if (gff) {
        def ch_gff = gff.endsWith('.gz') 
            ? GUNZIP_GTF([ [:], file(gff, checkIfExists: true) ]).gunzip
            : Channel.value([ [:], file(gff, checkIfExists: true) ])
        ch_gtf = GFFREAD(ch_gff, []).gtf.map { it[1] }
    }

    // Handle FASTA
    def fasta_provided = fasta ? true : false
    ch_fasta = fasta_provided 
        ? (fasta.endsWith('.gz') 
            ? GUNZIP_FASTA([ [:], file(fasta, checkIfExists: true) ]).gunzip.map { it[1] } 
            : Channel.value(file(fasta, checkIfExists: true)))
        : Channel.empty()

    // Generate chromosome sizes from FASTA (needed for Kallisto)
    ch_chrom_sizes = Channel.empty()
    if (fasta_provided) {
        ch_fasta_meta = ch_fasta.map { [ [:], it ] }
        CUSTOM_GETCHROMSIZES(ch_fasta_meta)
        ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
        ch_versions = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)
    }

    // Filter GTF if needed
    def filter_gtf_needed = (
        pseudo_aligner ||
        (!transcript_fasta)
    ) && !skip_gtf_filter

    if (filter_gtf_needed) {
        GTF_FILTER(ch_fasta, ch_gtf)
        ch_gtf = GTF_FILTER.out.genome_gtf.first()
    }

    // Handle transcript FASTA
    ch_transcript_fasta = Channel.empty()
    if (transcript_fasta) {
        ch_transcript_fasta = transcript_fasta.endsWith('.gz') 
            ? GUNZIP_TRANSCRIPT_FASTA([ [:], file(transcript_fasta, checkIfExists: true) ]).gunzip.map { it[1] } 
            : Channel.value(file(transcript_fasta, checkIfExists: true))
        if (gencode) {
            PREPROCESS_TRANSCRIPTS_FASTA_GENCODE(ch_transcript_fasta)
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.fasta
        }
    } else if (fasta_provided) {
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA(ch_fasta, ch_gtf).transcript_fasta
    }

    // Kallisto index handling
    ch_index = Channel.empty()

    def use_kallisto = pseudo_aligner == 'kallisto'
    def index_name = index != null ? (index instanceof String ? index : index.name) : null
    def index_provided = index_name != null && (index_name.endsWith('.idx') || index_name.endsWith('.tar.gz'))

    if (index_provided) {
        if (index_name.endsWith('.tar.gz')) {
            def ch_untar = UNTAR_KALLISTO_INDEX([ [:], index ]).untar
            ch_index = ch_untar.map { [ [:], it ] }
            ch_versions = ch_versions.mix(UNTAR_KALLISTO_INDEX.out.versions)
        } else {
            ch_index = Channel.of([ [:], file(index) ])
        }
    } else if (use_kallisto) {
        ch_index = KALLISTO_INDEX(ch_transcript_fasta.map { [ [:], it ] }).index
        ch_versions = ch_versions.mix(KALLISTO_INDEX.out.versions)
    }

    emit:
    fasta            = ch_fasta
    gtf              = ch_gtf
    transcript_fasta = ch_transcript_fasta
    index            = ch_index
    chrom_sizes      = ch_chrom_sizes
    versions         = ch_versions.ifEmpty(null)
}
