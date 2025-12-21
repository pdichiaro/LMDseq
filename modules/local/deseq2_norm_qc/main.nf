/*
 * DESeq2 normalization and quality control analysis (all genes)
 * Adapted from nf-core/rnaseq NORMALIZE_DESEQ2_QC_ALL_GENES module
 */
process DESEQ2_NORM_QC {
    tag "deseq2_qc_all_genes"
    label "process_medium"

    container 'docker://pdichiaro/lmdseq:latest'

    input:
    path counts
    val quantifier

    output:
    path "scaling_dat.txt"       , emit: scaling_factors
    path "scaling_factors/*_scaling_factor.txt", emit: scaling_factors_individual
    path "*_normalized_counts.txt", emit: normalized_counts
    path "*_rlog_counts.txt"     , optional:true, emit: rlog_counts
    path "*.RData"               , optional:true, emit: rdata
    path "*.pca.vals.txt"        , optional:true, emit: pca_all_genes_txt
    path "*.pca.top*.vals.txt"   , optional:true, emit: pca_top_genes_txt
    path "*.sample.dists.txt"    , optional:true, emit: sample_distances_txt
    // sample.correlations.txt - not created by R script (correlation available in PDF/PNG only)
    path "*.read.distribution.normalized.txt", optional:true, emit: read_dist_norm_txt
    // read.distribution.raw.txt - not created by R script (only normalized version is saved)
    // histogram/density txt files - not created by R script (data embedded in plots only)
    path "Quality_Control"       , optional:true, emit: quality_control_plots
    path "Read_Distribution"     , optional:true, emit: read_dist_plots
    path "*.log"                 , optional:true, emit: log
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    def label_lower = args2.toLowerCase()
    def label_upper = args2.toUpperCase()
    prefix = task.ext.prefix ?: "deseq2_all_genes"
    """
    echo "=== DESEQ2 QC DEBUG INFO ==="
    echo "Count file input: $counts"
    echo "Working directory: \$(pwd)"
    echo "Files in directory:"
    ls -la
    echo "Count file exists: \$(test -f $counts && echo YES || echo NO)"
    echo "Count file size: \$(wc -l $counts 2>/dev/null || echo 'Cannot read')"
    echo "=========================="
    
    normalize_deseq2_qc_all_genes.r \\
        --count_file $counts \\
        --outdir ./ \\
        --outprefix "$prefix" \\
        --quantifier "$quantifier" \\
        --cores $task.cpus \\
        $args

    echo "=== FILES GENERATED FOR MULTIQC ==="
    ls -la *.txt 2>/dev/null || echo "No .txt files found"
    echo "===================================="

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