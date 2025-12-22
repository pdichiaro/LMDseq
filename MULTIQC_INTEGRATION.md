# MultiQC Integration for DESeq2 Script

## Overview
The `deseq2_norm_qc.R` script has been updated to generate MultiQC-compatible TSV files with "DESeq2" namespace for integration into MultiQC reports.

## Modifications Made

### 1. Added MultiQC Export Function
**Location:** After directory creation (line ~83)

```r
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
```

### 2. PCA Data Export
**Location:** After PCA plot generation (line ~440)

**Output File:** `8_Quality_folder/deseq2_pca_mqc.tsv`

**Content:**
- PC1 coordinates
- PC2 coordinates
- Biological replicates
- Includes variance percentages in description

**Code:**
```r
pca_export <- pcaData[, c("PC1", "PC2", "Bio_replicates")]
rownames(pca_export) <- pcaData$name
write_multiqc_tsv(
  pca_export,
  file.path(Quality_folder, "deseq2_pca_mqc.tsv"),
  "DESeq2 PCA",
  paste0("Principal Component Analysis - PC1: ", percentVar[1], "%, PC2: ", percentVar[2], "%")
)
```

### 3. Sample Distance Matrix Export
**Location:** After sample distance heatmap generation (line ~390)

**Output File:** `8_Quality_folder/deseq2_sample_distances_mqc.tsv`

**Content:**
- Sample-to-sample Euclidean distances
- Based on VST-transformed data
- Full distance matrix with sample names

**Code:**
```r
write_multiqc_tsv(
  sampleDistMatrix,
  file.path(Quality_folder, "deseq2_sample_distances_mqc.tsv"),
  "DESeq2 Sample Distances",
  "Sample-to-sample Euclidean distances (VST-transformed)"
)
```

### 4. Read Distribution Statistics Export
**Location:** After read distribution plots (line ~348)

**Output File:** `8_Quality_folder/deseq2_read_distribution_mqc.tsv`

**Content:**
- Total reads per sample
- Mean counts per sample
- Median counts per sample
- Number of genes detected per sample

**Code:**
```r
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
```

## Output Files Generated

All files are created in the `8_Quality_folder` directory with the suffix `_mqc.tsv`:

1. **deseq2_pca_mqc.tsv** - PCA coordinates and variance
2. **deseq2_sample_distances_mqc.tsv** - Sample distance matrix
3. **deseq2_read_distribution_mqc.tsv** - Read distribution statistics

## MultiQC Configuration

All generated files use:
- **Namespace:** `DESeq2`
- **Plot Type:** `table`
- **Custom ID:** `deseq2_custom_data`

## File Format

Each TSV file contains:
```
# plot_type: 'table'
# section_name: '<Section Name>'
# description: '<Description>'
# pconfig:
#     id: 'deseq2_custom_data'
#     namespace: 'DESeq2'
<TAB-DELIMITED DATA>
```

## Usage

The script automatically generates these files during normal execution:

```bash
Rscript deseq2_norm_qc.R \
  -i EX_reads_RAW.txt \
  -m Sample.txt \
  -o output_directory \
  --min_reads 10
```

## Integration with MultiQC

MultiQC will automatically detect and parse these files when running:

```bash
multiqc output_directory/
```

The data will appear in the MultiQC report under the "DESeq2" section.

## Notes

- All statistics are based on normalized and filtered data
- PCA uses VST-transformed counts
- Sample distances use Euclidean distance on VST-transformed data
- Read distribution statistics calculated on DESeq2-normalized counts
- Script automatically handles both single and multiple replicate scenarios (see DESEQ2_DESIGN_FIX.md)
