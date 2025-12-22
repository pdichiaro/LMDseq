# DESeq2 Design Formula Strategy

## Background

The DESeq2 normalization and QC script uses a fixed `design ~ 1` formula, following the same approach as `nf-core/rnaseq` and the `pdichiaro/rnaseq` pipeline.

## Why `design ~ 1`?

Using an intercept-only model (`design ~ 1`) for normalization and quality control provides several advantages:

1. **Universal compatibility**: Works with any sample configuration (single or multiple replicates)
2. **Blind transformation**: Equivalent to setting `blind=TRUE` for variance-stabilizing transformation
3. **Focus on QC**: Appropriate for normalization and exploratory analysis without assuming experimental design
4. **Avoids errors**: Prevents issues when samples have identical replicate values
5. **Consistency**: Matches the strategy used in nf-core/rnaseq pipeline

## Previous Problem (Now Fixed)

Previously, the script used `design ~ Bio_replicates`, which caused an error when all samples had the same replicate value:

```
Error in DESeqDataSet(se, design = design, ignoreRank) : 
  design has a single variable, with all samples having the same value.
  use instead a design of '~ 1'.
```

## Current Implementation

The script now **always** uses `design ~ 1` for normalization and QC purposes.

### Code Changes

**Location:** `bin/deseq2_norm_qc.R` (line ~230)

```r
# Create DESeq2 dataset with design ~ 1
# Using intercept-only model (design=~1) for normalization and QC, consistent with nf-core/rnaseq approach.
# This is equivalent to setting blind=TRUE for variance-stabilizing transformation.
# This design works for all sample configurations (single or multiple replicates) and is appropriate
# for normalization, size factor estimation, and quality control analyses.
cat("Creating DESeq2 dataset with design ~ 1 for normalization and QC\n")
dds_ex_test <- DESeqDataSetFromMatrix(countData = matrix_test, colData = colData_test, design = ~1)
```

## Behavior

The script **always** uses `design ~ 1` regardless of the sample configuration.

### Console Output

```
Creating DESeq2 dataset with design ~ 1 for normalization and QC
```

### Results

- ✅ Works with single replicate samples
- ✅ Works with multiple replicate samples
- ✅ Performs blind variance-stabilizing transformation
- ✅ No errors regardless of replicate structure

## Impact

The `design ~ 1` approach is ideal for normalization and quality control purposes:

### What Works Perfectly

All normalization and QC functions work with `design ~ 1`:
- ✅ Size factor estimation (`estimateSizeFactors`)
- ✅ Variance stabilizing transformation (`vst`) with blind=TRUE behavior
- ✅ rlog transformation
- ✅ PCA analysis
- ✅ Sample distance calculations
- ✅ Read distribution statistics
- ✅ MultiQC export files
- ✅ Exploratory data analysis
- ✅ Works with any sample configuration

### When to Use Different Designs

- **`design ~ 1`**: For normalization, QC, and exploratory analysis (current implementation)
- **`design ~ condition`**: For differential expression analysis (separate DESeq2 workflow)

### Note on Differential Expression

The normalization script focuses on QC and data normalization. For differential expression analysis:
- Use a separate DESeq2 workflow with appropriate design formula
- The normalized counts from this script can be used as input
- Design formula should reflect the experimental conditions to be compared

## Usage

No changes required in pipeline usage. The script automatically detects and handles both cases:

```bash
# Works for both single and multiple replicate scenarios
nextflow run main.nf --input samplesheet.csv --outdir results
```

## Testing

The fix has been validated to:
1. ✅ Parse correctly in R
2. ✅ Handle single replicate samples without error
3. ✅ Maintain backward compatibility with multiple replicate samples
4. ✅ Generate all expected output files in both cases

## Reference Implementation

This approach is consistent with:
- **nf-core/rnaseq**: Uses `design ~ 1` for normalization modules
- **pdichiaro/rnaseq**: Both `normalize_deseq2_qc_all_genes.r` and `normalize_deseq2_qc_invariant_genes.r` use `design ~ 1`

From `pdichiaro/rnaseq/bin/normalize_deseq2_qc_all_genes.r`:
```r
# `design=~1` creates intercept-only model, equivalent to setting `blind=TRUE` for transformation.
dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=coldata, design=~1)
```

## Commit Information

**Initial Fix:** `1c7fa97` (conditional design based on replicates)  
**Updated to design ~ 1:** TBD (following nf-core/rnaseq strategy)  
**Date:** 2025-12-22  
**Branch:** main
