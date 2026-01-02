#!/usr/bin/env Rscript

#
# Created by Pierluigi Di Chiaro
# DESeq2 Normalization and Quality Control Script
# Licensed under MIT License - see LICENSE file for details
#

library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input gene expression matrix (EX_reads_RAW.txt)", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Sample metadata file (Sample.txt)", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="./",
              help="Output directory", metavar="character"),
  make_option(c("--min_reads"), type="integer", default=10,
              help="Minimum read count threshold [default %default]", metavar="integer")
)

optparser = OptionParser(option_list=option_list);
opt = parse_args(optparser)

if (is.null(opt$input)) {
  print_help(optparser)
  stop("Input gene expression matrix must be supplied", call.=FALSE)
}
if (is.null(opt$metadata)) {
  print_help(optparser)
  stop("Sample metadata file must be supplied", call.=FALSE)
}

# Load required libraries
suppressPackageStartupMessages({
  library("tidyverse")
  library("Rfast")
  library("DESeq2")
  library("Biobase")
  library("RColorBrewer")
  library("factoextra")
  library("gplots")
  library("ggrepel")
  library("pheatmap")
  library("PoiClaClu")

  library("fitdistrplus")
  library("FNN")
  library("igraph")
  library("dendsort")
  library("viridis")
  library("uwot")
})

# Set parameters
input_file <- opt$input
metadata_file <- opt$metadata
output_dir <- opt$outdir
min_reads <- opt$min_reads

# Hardcoded gene types to filter
gene_types <- c("protein_coding", "protein-coding", "protein-coding gene")

cat("Starting DESeq2 normalization and QC analysis...\n")
cat("Input file:", input_file, "\n")
cat("Metadata file:", metadata_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Minimum reads threshold:", min_reads, "\n")
cat("Gene types to keep:", paste(gene_types, collapse=", "), "\n")

# Create output directories
Norm_folder <- file.path(output_dir, "6_Norm_folder")
Counts_folder <- file.path(output_dir, "7_Counts_folder")
Quality_folder <- file.path(output_dir, "8_Quality_folder")

dir.create(Norm_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(Counts_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(Quality_folder, showWarnings = FALSE, recursive = TRUE)

# Function to write MultiQC TSV files with DESeq2 header
write_multiqc_tsv <- function(data, filename, section_name, description) {
  # Create header lines
  header_lines <- c(
    "# plot_type: 'table'",
    paste0("# section_name: '", section_name, "'"),
    paste0("# description: '", description, "'"),
    "# pconfig:",
    "#     id: 'deseq2_custom_data'",
    "#     namespace: 'DESeq2'"
  )
  
  # Write header
  writeLines(header_lines, filename)
  
  # Append data table
  write.table(data, filename, append = TRUE, sep = "\t", 
              quote = FALSE, row.names = TRUE, col.names = NA)
}

# Read metadata file (Sample.txt)
cat("Reading metadata file...\n")
MASTER_FILE <- read.delim(metadata_file, header=TRUE)

# Process sample names
MASTER_FILE$SHORT_NAME <- gsub(" ","", paste0(MASTER_FILE$SHORT_NAME, "_", MASTER_FILE$REPLICATE))
rownames(MASTER_FILE) <- MASTER_FILE$SHORT_NAME 
samples <- as.character(MASTER_FILE$SHORT_NAME)

# Create colData for DESeq2 - only sample_name and Bio_replicates
Bio_replicates <- as.character(MASTER_FILE$REPLICATE)

colData <- data.frame(
  sample_name = samples,
  Bio_replicates = Bio_replicates,
  row.names = samples
)

cat("Sample information:\n")
print(colData)

# Read gene expression matrix
cat("Reading gene expression matrix:", input_file, "\n")
all.reads <- read.delim(file=input_file, sep="\t", row.names=1, check.names=FALSE)
colnames(all.reads) <- gsub("^X","", colnames(all.reads))

# Check if this is minimal format (gene_id only) or full annotation format
# Full annotation format has annotation columns like "Symbol", "type_of_gene", etc.
# Minimal format only has sample columns (gene_id is in rownames)
has_annotations <- any(c("Symbol", "type_of_gene", "seqnames", "start") %in% colnames(all.reads))
cat("Matrix columns:", paste(colnames(all.reads), collapse=", "), "\n")
cat("Expected samples:", paste(samples, collapse=", "), "\n")
cat("Has annotations:", has_annotations, "\n")

if (has_annotations) {
  cat("Full annotation format detected\n")
  # Extract annotation columns (first 15 columns contain gene information)
  all_name <- all.reads[,1:15]
  count_data <- all.reads[, samples, drop=FALSE]
  count_data <- round(count_data)
  # Combine annotation with counts
  all.reads <- cbind(all_name, count_data)
} else {
  cat("Minimal format detected - gene_id only\n")
  # Create minimal annotation structure
  gene_ids <- rownames(all.reads)
  count_data <- all.reads[, samples, drop=FALSE]
  count_data <- round(count_data)
  
  # Create basic annotation columns for compatibility
  all_name <- data.frame(
    gene_ids = gene_ids,
    IDs = gene_ids,
    tx_id = seq_along(gene_ids),
    seqnames = "unknown",
    start = "unknown",
    end = "unknown", 
    width = "unknown",
    strand = "unknown",
    LOC = "unknown",
    type_of_gene = "protein_coding",  # Default to protein_coding for filtering
    Symbol = gene_ids,  # Use gene_id as symbol
    HGNC = "unknown",
    Name = "unknown", 
    Entrez.Gene.ID = "unknown",
    cDNA_LENGTH = "unknown"
  )
  
  # Combine annotation with counts  
  all.reads <- cbind(all_name, count_data)
}

# Filter by gene type
cat("Filtering genes by type...\n")
cat("Available gene types:", unique(all.reads$type_of_gene), "\n")

ww <- which(all.reads$type_of_gene %in% gene_types)
if(length(ww) == 0) {
  cat("Warning: No genes found with specified types. Using all genes.\n")
  ww <- 1:nrow(all.reads)
}
all.reads <- all.reads[ww,]

cat("Genes after type filtering:", nrow(all.reads), "\n")

# Filter low-expressed genes
cat("Filtering low-expressed genes...\n")
all.reads_c <- all.reads

# Use Symbol column as rownames (duplicates should be handled upstream in create_reference_db.R)
cat("Setting rownames using Symbol column...\n")

# Safety check for remaining duplicates
symbols_to_use <- all.reads$Symbol
if(any(duplicated(symbols_to_use))) {
  cat("Warning: Found", sum(duplicated(symbols_to_use)), "duplicate symbols - this should have been handled upstream!\n")
  cat("Duplicated symbols:", paste(unique(symbols_to_use[duplicated(symbols_to_use)]), collapse=", "), "\n")
  stop("Duplicate symbols detected - please check create_reference_db.R output")
} else {
  cat("No duplicate symbols found - proceeding with Symbol rownames\n")
}

rownames(all.reads_c) <- symbols_to_use

# Create binary matrix for filtering
Rep_counts_all <- as.data.frame(matrix(0, nrow=nrow(all.reads_c), ncol=length(samples)))
colnames(Rep_counts_all) <- samples
rownames(Rep_counts_all) <- rownames(all.reads_c)

for(c in colnames(Rep_counts_all)){
  w <- which(all.reads_c[,c] >= min_reads)
  Rep_counts_all[w,c] <- 1
}

# Filter genes that have sufficient expression in at least half of samples
keep_genes <- rowSums(Rep_counts_all) >= round(length(samples)/2)
all.reads_c <- all.reads_c[keep_genes,]
rownames(all.reads_c) <- NULL
all_name <- all.reads_c[,1:15]

cat("Genes after expression filtering:", nrow(all.reads_c), "\n")

# Save filtered counts
write.table(all.reads_c, file=file.path(Counts_folder, "EX_reads_RAW_filt.txt"), sep="\t", col.names=NA)

# Prepare data for DESeq2
cat("Running DESeq2 normalization...\n")
RAED_MAT <- all.reads_c
colData_test <- colData[samples,, drop=FALSE]
matrix_test <- as.matrix(RAED_MAT[,rownames(colData_test)])

# Note: Integer conversion already handled upstream in create_reference_db.R

# Create DESeq2 dataset with design ~ 1
# Using intercept-only model (design=~1) for normalization and QC, consistent with nf-core/rnaseq approach.
# This is equivalent to setting blind=TRUE for variance-stabilizing transformation.
# This design works for all sample configurations (single or multiple replicates) and is appropriate
# for normalization, size factor estimation, and quality control analyses.
cat("Creating DESeq2 dataset with design ~ 1 for normalization and QC\n")
dds_ex_test <- DESeqDataSetFromMatrix(countData = matrix_test, colData = colData_test, design = ~1)
dds_ex_test <- estimateSizeFactors(dds_ex_test)

# Save normalization parameters
size <- as.data.frame(sizeFactors(dds_ex_test))                                      
colnames(size) <- "Depth"
write.table(size, file=file.path(Counts_folder, "Normalisation_Parameters.txt"), sep="\t", col.names=NA)

# Save scaling factors in nf-core format
scaling_dat <- data.frame(
  sample = names(sizeFactors(dds_ex_test)),
  size_factor = as.numeric(sizeFactors(dds_ex_test))
)
write.table(scaling_dat, file=file.path(output_dir, "scaling_dat.txt"), 
            sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# Create scaling_factors directory and save individual files
scaling_factors_dir <- file.path(output_dir, "scaling_factors")
dir.create(scaling_factors_dir, showWarnings = FALSE, recursive = TRUE)
for (samplename in names(sizeFactors(dds_ex_test))) {
  factor_file <- file.path(scaling_factors_dir, paste0(samplename, "_scaling_factor.txt"))
  write(as.numeric(sizeFactors(dds_ex_test)[samplename]), file=factor_file)
}

# Get normalized counts
all.reads_d <- as.data.frame(counts(dds_ex_test, normalized = TRUE)) 
all.reads_d <- cbind(all_name, all.reads_d)

# Save normalized counts
write.table(all.reads_d, file=file.path(Counts_folder, "EX_reads_NORM_filt.txt"), sep="\t", col.names=NA)

# Plot read distribution - Raw data
cat("Creating read distribution plots...\n")
all.reads_t <- as.data.frame(log(all.reads[,samples]+1))
all.reads_tt <- all.reads_t[,rownames(colData)] %>% rownames_to_column(var="gene") %>% as_tibble()            
gathered_all.reads <- gather(all.reads_tt, key = "samplename", value = "normalized_counts", -gene)
gathered_all.reads <- gathered_all.reads %>% 
  left_join(colData %>% rownames_to_column("samplename") %>% dplyr::select(samplename, Bio_replicates), by = c("samplename")) %>%
  arrange(Bio_replicates) %>% 
  mutate(samplename = factor(samplename, levels = unique(samplename)))

nm <- file.path(Norm_folder, "Read_Distribution_Raw.pdf")
pdf(nm, width=20, height=20)
par(mfrow=c(1,3))

p.1 <- ggplot(gathered_all.reads, aes(x = samplename, y = normalized_counts, fill = Bio_replicates)) + 
  geom_boxplot(show.legend = FALSE) + xlab("") +
  ylab(expression(ln(count + 1))) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.2 <- ggplot(gathered_all.reads, aes(x = normalized_counts, colour = Bio_replicates, fill = Bio_replicates)) +
  geom_histogram(binwidth = 1) + xlab(expression(ln(count + 1))) + ylab("frequency") + 
  ylim(c(0,200000)) + theme(legend.position = "top") + theme_classic()

p.3 <- ggplot(gathered_all.reads, aes(x = normalized_counts, colour = Bio_replicates, fill = Bio_replicates)) + 
  geom_density(alpha = 0.2, size = 1.25) + xlab(expression(ln(count))) + ylim(c(0, 0.5)) +
  theme(legend.position = "top") + theme_classic()

print(p.1)
print(p.2)
print(p.3)
dev.off()

# Plot read distribution - Normalized data
all.reads_z <- as.data.frame(log(all.reads_d[,samples]+1))
all.reads_zz <- all.reads_z[,rownames(colData)] %>% rownames_to_column(var="gene") %>% as_tibble()            
gathered_all.reads <- gather(all.reads_zz, key = "samplename", value = "normalized_counts", -gene)
gathered_all.reads <- gathered_all.reads %>% 
  left_join(colData %>% rownames_to_column("samplename") %>% dplyr::select(samplename, Bio_replicates), by = c("samplename")) %>%
  arrange(Bio_replicates) %>% 
  mutate(samplename = factor(samplename, levels = unique(samplename)))

nm <- file.path(Norm_folder, "Read_Distribution_Norm_Filt.pdf")
pdf(nm, width=10, height=10)
par(mfrow=c(1,3))

p.1 <- ggplot(gathered_all.reads, aes(x = samplename, y = normalized_counts, fill = Bio_replicates)) + 
  geom_boxplot(show.legend = FALSE) + xlab("") +
  ylab(expression(ln(count + 1))) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.2 <- ggplot(gathered_all.reads, aes(x = normalized_counts, colour = Bio_replicates, fill = Bio_replicates)) +
  geom_histogram(binwidth = 1) + xlab(expression(ln(count + 1))) + ylab("frequency") + 
  ylim(c(0,200000)) + theme(legend.position = "top") + theme_classic()

p.3 <- ggplot(gathered_all.reads, aes(x = normalized_counts, colour = Bio_replicates, fill = Bio_replicates)) + 
  geom_density(alpha = 0.2, size = 1.25) + xlab(expression(ln(count))) + ylim(c(0, 0.5)) +
  theme(legend.position = "top") + theme_classic()

print(p.1)
print(p.2)
print(p.3)
dev.off()

# Export read distribution statistics for MultiQC
read_stats <- data.frame(
  Total_Reads = colSums(all.reads_d[, samples]),
  Mean_Counts = colMeans(all.reads_d[, samples]),
  Median_Counts = apply(all.reads_d[, samples], 2, median),
  Genes_Detected = colSums(all.reads_d[, samples] > 0)
)
write_multiqc_tsv(
  read_stats,
  file.path(Quality_folder, "deseq2_read_distribution_mqc.tsv"),
  "DESeq2 Read Distribution",
  "Read distribution statistics for normalized and filtered data"
)

# Quality control analysis
cat("Performing quality control analysis...\n")
dir.create(Quality_folder, showWarnings = FALSE, recursive = TRUE)

# Variance stabilizing transformation
vsd <- vst(dds_ex_test, blind = TRUE)

# Save rlog transformed data
rlog <- cbind(all_name, as.data.frame(assay(vsd)))
write.table(rlog, file=file.path(Counts_folder, "rLog_reads.txt"), sep="\t", col.names=NA)

# Create annotation for heatmaps
annot_c <- as.data.frame(matrix(NA, nrow=length(rownames(colData)), ncol=1))
rownames(annot_c) <- rownames(colData)
colnames(annot_c) <- c("Bio_replicates")
annot_c$Bio_replicates <- colData$Bio_replicates

# Sample-to-sample distances
cat("Creating sample distance heatmaps...\n")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists, labels=TRUE)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pdf(file.path(Quality_folder, "Heatmap_sampleTosample_distances.pdf"))
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         fontsize_row=6,
         fontsize_col=6,
         border_color = FALSE,
         col = colors)
dev.off()

# Export sample distances for MultiQC
write_multiqc_tsv(
  sampleDistMatrix,
  file.path(Quality_folder, "deseq2_sample_distances_mqc.tsv"),
  "DESeq2 Sample Distances",
  "Sample-to-sample Euclidean distances (VST-transformed)"
)

# Poisson distance
poisd <- PoissonDistance(t(counts(dds_ex_test)))
samplePoisDistMatrix <- as.matrix(poisd$dd, labels=TRUE)
rownames(samplePoisDistMatrix) <- colnames(counts(dds_ex_test))[as.numeric(rownames(samplePoisDistMatrix))]
colnames(samplePoisDistMatrix) <- NULL

pdf(file.path(Quality_folder, "Heatmap_Poisson_Distance.pdf"))
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         fontsize_row=6,
         fontsize_col=6,
         border_color = FALSE,
         col = colors)
dev.off()

# Correlation heatmap
Pearson.corr <- cor(assay(vsd), method="pearson")  
pdf(file.path(Quality_folder, "Heatmap_sampleTosample_correlation.pdf"), height=10, width=10)
pheatmap(Pearson.corr,
         annotation_col = annot_c,
         fontsize = 6,
         fontsize_row=6,
         fontsize_col=6,
         border_color = FALSE,
         col=rev(brewer.pal(8, "RdBu")))                                      
dev.off()

# PCA analysis
cat("Performing PCA analysis...\n")
# Use all genes for PCA
ntop_pca <- nrow(matrix_test)
cat("Using all", ntop_pca, "genes for PCA\n")
pcaData <- plotPCA(vsd, intgroup = c("Bio_replicates"), ntop = ntop_pca, returnData=TRUE)
percentVar <- round(100*attr(pcaData,"percentVar"))

# Create PCA plot with proper dimensions
gg <- ggplot(pcaData, aes(PC1, PC2, color=Bio_replicates, label=name)) +
  geom_point(size=3) +
  geom_text_repel(max.overlaps = Inf, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() +
  theme(aspect.ratio = 1)

# Save with explicit dimensions to avoid viewport errors
ggsave("PCA_rlogTransformed.pdf", plot = gg, device="pdf", useDingbats=FALSE, 
       width=10, height=10, path=Quality_folder, limitsize = FALSE, units = "in")

# Export PCA data for MultiQC
pca_export <- pcaData[, c("PC1", "PC2", "Bio_replicates")]
rownames(pca_export) <- pcaData$name
write_multiqc_tsv(
  pca_export,
  file.path(Quality_folder, "deseq2_pca_mqc.tsv"),
  "DESeq2 PCA",
  paste0("Principal Component Analysis - PC1: ", percentVar[1], "%, PC2: ", percentVar[2], "%")
)

# Extended PCA analysis
# Use all genes (consistent with main PCA analysis)
ntop = nrow(assay(vsd))
intgroup = c("Bio_replicates")
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
centered_scaled <- scale(t(assay(vsd[select,])), center=TRUE, scale=FALSE)
pca <- prcomp(centered_scaled, scale. = FALSE) 
percentVar <- round(100*(pca$sdev^2/sum(pca$sdev^2)))
intgroup.df <- as.data.frame(colData(vsd)[, intgroup, drop = FALSE])
d <- data.frame(pca$x, intgroup.df, name = colnames(vsd))
attr(d, "percentVar") <- percentVar[1:2]

# Percent_var_PCA, 3D PCA, and corrplot removed per user request 

cat("DESeq2 normalization and QC analysis completed successfully!\n")
cat("Output files saved to:", output_dir, "\n")
