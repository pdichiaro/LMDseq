/*
 * DESeq2 normalization and quality control analysis
 */
process DESEQ2_NORM_QC {
    tag "deseq2_norm_qc"
    label "process_medium"

    container 'docker://pdichiaro/lmdseq:latest'

    input:
    path gene_matrix
    path sample_metadata

    output:
    path "6_Norm_folder/**"                        , emit: norm_files
    path "7_Counts_folder/**"                      , emit: count_files
    path "8_Quality_folder/**"                     , emit: qc_files
    path "deseq2_normalized_counts.txt", optional: true, emit: normalized_counts
    path "deseq2_rlog_counts.txt", optional: true , emit: rlog_counts
    path "deseq2.dds.RData", optional: true        , emit: dds_rdata
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "=== DESEQ2 QC DEBUG INFO ==="
    echo "Count file input: $gene_matrix"
    echo "Metadata file input: $sample_metadata"
    echo "Working directory: \$(pwd)"
    echo "Files in directory:"
    ls -la
    echo "Count file exists: \$(test -f $gene_matrix && echo YES || echo NO)"
    echo "Metadata file exists: \$(test -f $sample_metadata && echo YES || echo NO)"
    echo "=========================="
    
    Rscript ${projectDir}/bin/deseq2_norm_qc.R \\
        --input $gene_matrix \\
        --metadata $sample_metadata \\
        --outdir . \\
        --min_reads ${params.min_reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """

    stub:
    def args2 = task.ext.args2 ?: ''
    def label_lower = args2.toLowerCase()
    prefix = task.ext.prefix ?: "deseq2_all_genes"
    """
    touch ${prefix}.dds.RData
    touch ${prefix}_rlog_counts.txt
    touch ${prefix}.pca.vals.txt
    touch ${prefix}.pca.top500.vals.txt
    touch ${prefix}.sample.dists.txt
    touch ${prefix}.read.distribution.normalized.txt
    touch ${prefix}_normalized_counts.txt
    touch scaling_dat.txt
    touch R_sessionInfo.log
    
    mkdir -p Quality_Control
    touch Quality_Control/Sample_Distance_Heatmap.pdf
    touch Quality_Control/Sample_Distance_Heatmap.png
    touch Quality_Control/Sample_Correlation_Heatmap.pdf
    touch Quality_Control/Sample_Correlation_Heatmap.png
    touch Quality_Control/PCA_Plot_All_Genes.pdf
    touch Quality_Control/PCA_Plot_Top500_Genes.pdf
    
    mkdir -p Read_Distribution
    touch Read_Distribution/Read_Distribution_Raw_boxplot.pdf
    touch Read_Distribution/Read_Distribution_Raw_histogram.pdf
    touch Read_Distribution/Read_Distribution_Raw_density.pdf
    touch Read_Distribution/Read_Distribution_Norm_Filt_boxplot.pdf
    touch Read_Distribution/Read_Distribution_Norm_Filt_histogram.pdf
    touch Read_Distribution/Read_Distribution_Norm_Filt_density.pdf

    mkdir -p scaling_factors
    touch scaling_factors/sample1_scaling_factor.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}