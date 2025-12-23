# PCA 3D and Corrplot Error Fix

**Date:** 2025-12-23  
**Commit:** `f79e270`  
**Issue:** Errors when analyzing datasets with few samples (< 3 samples)

---

## Problem Description

The DESeq2 script was failing when generating 3D PCA plots and correlation plots with the following errors:

### Error 1: 3D PCA Plot
```
Warning: Could not create 3D PCA plots: object 'PC3' not found
```

**Cause:** When there are fewer than 3 samples, PCA generates fewer than 3 principal components. The script was trying to access `PC3` which doesn't exist.

### Error 2: Corrplot
```
Error in r[i1] - r[-length(r):-(length(r) - lag + 1L)] : 
  non-numeric argument to binary operator
Calls: corrplot -> diff -> diff.default
Execution halted
```

**Cause:** The `corrplot()` function was receiving a dataframe containing non-numeric columns (`Bio_replicates`, `name`) mixed with PC columns, causing computation errors.

---

## Root Causes

### 1. No Check for PC3 Existence

**Previous code:**
```r
# 3D PCA plots
tryCatch({
  p <- plot_ly(d, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Bio_replicates, ...)
  # ...
}, error = function(e) {
  cat("Warning: Could not create 3D PCA plots:", e$message, "\n")
})
```

The `tryCatch` caught the error but **still caused pipeline failure** because the error occurred in critical plotting functions.

### 2. Non-Numeric Columns in Corrplot

**Previous code:**
```r
# PCA correlation plot
corrplot(as.matrix(d[,1:min(ncol(d), 10)]), ...)
```

The dataframe `d` structure:
```
d = [PC1, PC2, PC3?, ..., Bio_replicates, name]
     ^^^^^^^^^^^^^^^      ^^^^^^^^^^^^^^^^^^^^
     Numeric cols         Non-numeric cols
```

Selecting by position `d[,1:10]` included non-numeric columns, causing `corrplot()` to fail.

---

## Solution Implemented

### 1. Check PC3 Existence Before 3D Plot

```r
# 3D PCA plots (only if PC3 exists)
if("PC3" %in% colnames(d)) {
  tryCatch({
    p <- plot_ly(d, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Bio_replicates, colors = "Set1", 
                 marker = list(size=10, opacity=0.5))
    saveWidget(p, file=file.path(Quality_folder, "3D_PCA_Bio_replicates.html"), selfcontained=FALSE)
    
    p <- plot_ly(d, x = ~PC1, y = ~PC2, z = ~PC3, color = ~name, colors = "Set1", 
                 marker = list(size=10, opacity=0.5))
    saveWidget(p, file=file.path(Quality_folder, "3D_PCA_Samples.html"), selfcontained=FALSE)
  }, error = function(e) {
    cat("Warning: Could not create 3D PCA plots:", e$message, "\n")
  })
} else {
  cat("Skipping 3D PCA plots: Less than 3 principal components available\n")
}
```

**Key improvements:**
- ✅ Check if `PC3` exists before attempting 3D plot
- ✅ Graceful skip with informative message
- ✅ No pipeline failure

### 2. Select Only PC Columns for Corrplot

```r
# PCA correlation plot - use only PC columns (numeric)
pc_cols <- grep("^PC[0-9]+$", colnames(d), value = TRUE)
if(length(pc_cols) > 1) {
  pdf(file.path(Quality_folder, "All_pca.pdf"))
  corrplot(as.matrix(d[, pc_cols[1:min(length(pc_cols), 10)]]), is.corr=FALSE, 
           tl.cex = 0.7, cl.cex = 0.4, rect.lwd = 0.1, cl.pos = "b")
  dev.off()
} else {
  cat("Skipping PCA correlation plot: Not enough principal components\n")
}
```

**Key improvements:**
- ✅ Use `grep("^PC[0-9]+$", ...)` to select **only PC columns** (numeric)
- ✅ Exclude `Bio_replicates` and `name` columns
- ✅ Check that at least 2 PCs exist
- ✅ Graceful skip with informative message

---

## Benefits

### Robustness
- ✅ **No pipeline failures** with small sample sizes
- ✅ Handles edge cases (2 samples, 3 samples, etc.)
- ✅ Informative messages about what was skipped

### Correctness
- ✅ Only plots what's available
- ✅ Uses only numeric data for correlation plots
- ✅ No forced plotting of non-existent components

### User Experience
- ✅ Clear console messages
- ✅ Pipeline completes successfully
- ✅ Available plots are generated correctly

---

## Scenarios

| Samples | PCs Generated | 3D PCA | Corrplot |
|---------|---------------|--------|----------|
| 2 | 2 (PC1, PC2) | ⏭️ Skipped | ✅ Generated |
| 3 | 3 (PC1, PC2, PC3) | ✅ Generated | ✅ Generated |
| 4+ | 4+ | ✅ Generated | ✅ Generated |

---

## Console Output

### With < 3 Samples
```
Performing PCA analysis...
Using all 14193 genes for PCA
[PCA plot generated successfully]
Skipping 3D PCA plots: Less than 3 principal components available
[Corrplot generated successfully with available PCs]
```

### With ≥ 3 Samples
```
Performing PCA analysis...
Using all 14193 genes for PCA
[PCA plot generated successfully]
[3D PCA plots generated successfully]
[Corrplot generated successfully]
```

---

## Technical Details

### PC Column Detection
```r
# Regex pattern matches: PC1, PC2, PC3, ..., PC99, etc.
# Does NOT match: Bio_replicates, name, other metadata
pc_cols <- grep("^PC[0-9]+$", colnames(d), value = TRUE)
```

### Safe Column Selection
```r
# Before: d[,1:10] could include non-numeric columns
# After:  d[, pc_cols[1:min(length(pc_cols), 10)]]
#         ^^^^^^^^^^ Only numeric PC columns
```

---

## Related Files

- `bin/deseq2_norm_qc.R` - Main script (lines ~488-515)
- `PCA_FIX.md` - Main PCA viewport error fix
- `DESEQ2_DESIGN_FIX.md` - Design formula fix

---

## Testing

The fix has been validated to:
- ✅ Work with 2 samples (minimum)
- ✅ Work with 3+ samples (full PCA)
- ✅ Generate corrplot with only numeric data
- ✅ Complete pipeline without errors
- ✅ Provide clear console feedback

---

## Impact on Results

The fix:
- ✅ **Does NOT** change PCA computation
- ✅ **Does NOT** affect normalization or counts
- ✅ **Only affects** plot generation logic
- ✅ **Prevents** pipeline failures
- ✅ **Improves** robustness for small datasets

---

## Best Practices Applied

1. **Check before use**: Verify PC3 exists before plotting
2. **Select by name**: Use column names, not positions
3. **Type safety**: Ensure only numeric columns for numeric operations
4. **Graceful degradation**: Skip unavailable plots, continue pipeline
5. **User feedback**: Informative messages about skipped steps
