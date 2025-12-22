# DESeq2 Design Formula Fix

## Problem Description

When running the DESeq2 normalization script with samples that all have the same replicate value (e.g., all samples are "replicate 1"), DESeq2 throws an error:

```
Error in DESeqDataSet(se, design = design, ignoreRank) : 
  design has a single variable, with all samples having the same value.
  use instead a design of '~ 1'. estimateSizeFactors, rlog and the VST can then be used
```

This happens because DESeq2's design formula `~ Bio_replicates` cannot estimate parameters when there is no variation in the biological replicate variable.

## Solution Implemented

The script now automatically detects whether samples have multiple replicate values or a single replicate value, and adjusts the design formula accordingly.

### Code Changes

**Location:** `bin/deseq2_norm_qc.R` (line ~230)

```r
# Check if all samples have the same Bio_replicates value
unique_replicates <- unique(colData_test$Bio_replicates)
cat("Unique replicate values:", paste(unique_replicates, collapse=", "), "\n")

# Create DESeq2 dataset with appropriate design
if (length(unique_replicates) == 1) {
  cat("All samples have the same replicate value - using design ~ 1\n")
  dds_ex_test <- DESeqDataSetFromMatrix(countData = matrix_test, colData = colData_test, design = ~1)
} else {
  cat("Multiple replicate values detected - using design ~ Bio_replicates\n")
  dds_ex_test <- DESeqDataSetFromMatrix(countData = matrix_test, colData = colData_test, design = ~Bio_replicates)
}
```

## Behavior

### Case 1: Single Replicate Value (All samples identical)
**Example:** All samples have `Bio_replicates = "1"`

- **Design used:** `~ 1`
- **Console output:** 
  ```
  Unique replicate values: 1
  All samples have the same replicate value - using design ~ 1
  ```
- **Result:** Normalization proceeds successfully without error

### Case 2: Multiple Replicate Values
**Example:** Samples have `Bio_replicates = "1", "2", "3"`

- **Design used:** `~ Bio_replicates`
- **Console output:**
  ```
  Unique replicate values: 1, 2, 3
  Multiple replicate values detected - using design ~ Bio_replicates
  ```
- **Result:** Standard DESeq2 normalization with biological replicate information

## Impact

This fix allows the pipeline to handle both scenarios:

1. **Exploratory analysis** with single replicate samples (e.g., pilot studies, technical replicates only)
2. **Full differential expression analysis** with biological replicates

### What Still Works

All normalization and QC functions work with both designs:
- ✅ Size factor estimation (`estimateSizeFactors`)
- ✅ Variance stabilizing transformation (`vst`)
- ✅ rlog transformation
- ✅ PCA analysis
- ✅ Sample distance calculations
- ✅ Read distribution statistics
- ✅ MultiQC export files

### Limitations with `design ~ 1`

When using `design ~ 1` (single replicate case):
- ⚠️ Cannot perform differential expression analysis (no contrasts possible)
- ⚠️ Cannot estimate dispersion between biological replicates
- ✅ Can still perform normalization, quality control, and exploratory analysis

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

## Commit Information

**Commit Hash:** `1c7fa97`  
**Date:** 2025-12-22  
**Branch:** main
