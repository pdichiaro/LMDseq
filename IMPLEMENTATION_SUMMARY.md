# LMDseq Implementation Summary

**Last Updated:** 2025-12-23  
**Repository:** pdichiaro/LMDseq  
**Branch:** main

## Overview

This document summarizes the recent improvements made to the LMDseq pipeline, focusing on MultiQC integration, DESeq2 design formula strategy, and PCA optimization.

## Changes Implemented

### 1. MultiQC Integration ✅

**Commit:** `dac769a`

Added automatic generation of MultiQC-compatible TSV files for DESeq2 quality control metrics.

#### Features Added:
- **PCA Export**: `deseq2_pca_mqc.tsv`
  - PC1 and PC2 coordinates
  - Biological replicate information
  - Variance percentages in description

- **Sample Distance Matrix**: `deseq2_sample_distances_mqc.tsv`
  - Sample-to-sample Euclidean distances
  - Based on VST-transformed data

- **Read Distribution Statistics**: `deseq2_read_distribution_mqc.tsv`
  - Total reads per sample
  - Mean and median counts
  - Number of genes detected

#### Implementation:
```r
write_multiqc_tsv <- function(data, filename, section_name, description) {
  header_lines <- c(
    "# plot_type: 'table'",
    paste0("# section_name: '", section_name, "'"),
    paste0("# description: '", description, "'"),
    "# pconfig:",
    "#     id: 'deseq2_custom_data'",
    "#     namespace: 'DESeq2'"
  )
  writeLines(header_lines, filename)
  write.table(data, filename, append = TRUE, sep = "\t", 
              quote = FALSE, row.names = TRUE, col.names = NA)
}
```

**Documentation:** `MULTIQC_INTEGRATION.md`

---

### 2. DESeq2 Design Formula Strategy ✅

**Commits:** `1c7fa97`, `74aa911`, `9e8ab6b`

Changed DESeq2 design formula to use `design ~ 1` for all samples, following nf-core/rnaseq best practices.

#### Rationale:
1. **Universal Compatibility**: Works with any sample configuration
2. **Blind Transformation**: Equivalent to `blind=TRUE` for VST
3. **QC Focus**: Appropriate for normalization and exploratory analysis
4. **Consistency**: Matches nf-core/rnaseq and pdichiaro/rnaseq approach
5. **Error Prevention**: Avoids issues with identical replicate values

#### Previous Issue (Fixed):
```
Error: design has a single variable, with all samples having the same value.
use instead a design of '~ 1'.
```

#### Current Implementation:
```r
# Using intercept-only model (design=~1) for normalization and QC
cat("Creating DESeq2 dataset with design ~ 1 for normalization and QC\n")
dds_ex_test <- DESeqDataSetFromMatrix(countData = matrix_test, 
                                      colData = colData_test, 
                                      design = ~1)
```

#### Reference:
From `pdichiaro/rnaseq/bin/normalize_deseq2_qc_all_genes.r`:
```r
# `design=~1` creates intercept-only model, equivalent to blind=TRUE
dds <- DESeqDataSetFromMatrix(countData=round(counts), 
                               colData=coldata, 
                               design=~1)
```

**Documentation:** `DESEQ2_DESIGN_FIX.md`

---

### 3. PCA Optimization and Viewport Fix ✅

**Commit:** `472895e`, `1d3bff4`

Fixed viewport error during PCA plot generation and optimized PCA to use top 500 most variable genes.

#### Problem Solved:
```
Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : 
  Viewport has zero dimension(s)
```

Previously using all 14,193 genes which caused:
- Viewport calculation errors
- Slow computation
- Not following best practices

#### Changes Made:
- **Limit PCA to top 500 genes**: Standard practice, used by nf-core/rnaseq
- **Fix plot rendering**: Replace `coord_fixed()` with `theme(aspect.ratio = 1)`
- **Add max.overlaps**: Prevent label hiding in `geom_text_repel()`
- **Explicit dimensions**: Add `units = "in"` to `ggsave()`
- **Apply to both PCA sections**: Main and extended PCA analysis

#### Code:
```r
# Use top 500 most variable genes for PCA (consistent with nf-core/rnaseq)
ntop_pca <- min(500, nrow(matrix_test))
cat("Using top", ntop_pca, "most variable genes for PCA\n")
pcaData <- plotPCA(vsd, intgroup = c("Bio_replicates"), ntop = ntop_pca, returnData=TRUE)

# Create PCA plot with proper dimensions
gg <- ggplot(pcaData, aes(PC1, PC2, color=Bio_replicates, label=name)) +
  geom_point(size=3) +
  geom_text_repel(max.overlaps = Inf, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() +
  theme(aspect.ratio = 1)
```

#### Benefits:
- ✅ No more viewport errors
- ✅ Faster computation (500 vs 14,193 genes)
- ✅ Better biological interpretation
- ✅ Follows field standards
- ✅ Consistent with nf-core/rnaseq

**Documentation:** `PCA_FIX.md`

---

## Files Modified

### Core Scripts
- `bin/deseq2_norm_qc.R` - Main normalization and QC script

### Configuration
- `conf/modules.config` - Module configurations
- `modules/local/deseq2_norm_qc/main.nf` - Nextflow module
- `subworkflows/local/kallisto/main.nf` - Kallisto subworkflow
- `workflows/lmd_rnaseq.nf` - Main workflow

### Documentation (New)
- `MULTIQC_INTEGRATION.md` - MultiQC integration details
- `DESEQ2_DESIGN_FIX.md` - Design formula strategy explanation
- `PCA_FIX.md` - PCA optimization and viewport error fix
- `IMPLEMENTATION_SUMMARY.md` - This document

---

## Output Files Generated

All files are created in `8_Quality_folder/`:

### MultiQC TSV Files
1. `deseq2_pca_mqc.tsv`
2. `deseq2_sample_distances_mqc.tsv`
3. `deseq2_read_distribution_mqc.tsv`

### Standard Outputs (Unchanged)
- Normalized counts: `deseq2_normalized_counts.txt`
- rlog counts: `deseq2_rlog_counts.txt`
- Size factors: `scaling_factors/*.txt`
- DESeq2 object: `deseq2.dds.RData`
- Various QC plots (PDF)

---

## Compatibility

### Works With:
- ✅ Single replicate samples
- ✅ Multiple replicate samples
- ✅ Any experimental design
- ✅ Technical replicates only
- ✅ Biological replicates
- ✅ Mixed sample configurations

### MultiQC Integration:
- ✅ Automatic detection by MultiQC
- ✅ Displays in "DESeq2" section
- ✅ No additional configuration required

---

## Usage

No changes to pipeline usage required. Run as normal:

```bash
nextflow run pdichiaro/LMDseq \
  --input samplesheet.csv \
  --outdir results \
  --genome GRCh38
```

MultiQC will automatically include DESeq2 data:
```bash
multiqc results/
```

---

## Testing Status

- ✅ R syntax validation passed
- ✅ MultiQC TSV format validated
- ✅ Design ~ 1 tested with sample data
- ✅ Backward compatibility confirmed
- ✅ Documentation complete

---

## Benefits

1. **Robustness**: No more errors with identical replicates
2. **Consistency**: Matches nf-core standards
3. **Visibility**: DESeq2 QC metrics in MultiQC reports
4. **Flexibility**: Works with any sample configuration
5. **Documentation**: Comprehensive explanations provided

---

## Next Steps

The pipeline is now ready for production use with:
- Enhanced MultiQC reporting
- Robust DESeq2 normalization
- Universal sample compatibility

For differential expression analysis, use a separate workflow with appropriate design formula based on experimental conditions.

---

## References

- **nf-core/rnaseq**: https://nf-co.re/rnaseq
- **pdichiaro/rnaseq**: https://github.com/pdichiaro/rnaseq
- **DESeq2 documentation**: https://bioconductor.org/packages/DESeq2
- **MultiQC documentation**: https://multiqc.info/

---

## Contact

For questions or issues, please open an issue on the GitHub repository.
