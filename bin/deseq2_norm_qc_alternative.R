#!/usr/bin/env Rscript

#
# Created by Pierluigi Di Chiaro
# DESeq2 Normalization and Quality Control Script - Alternative version with unique symbols
# Licensed under MIT License - see LICENSE file for details
#

# This is an alternative fix that makes gene symbols unique instead of using gene_ids
# You can replace the original if you prefer to keep symbol-based rownames

# Function to make unique symbols
make_unique_symbols <- function(symbols, gene_ids) {
  # Find duplicated symbols
  duplicated_mask <- duplicated(symbols) | duplicated(symbols, fromLast = TRUE)
  
  if(sum(duplicated_mask) > 0) {
    cat("Making unique symbols for", sum(duplicated_mask), "duplicated entries...\n")
    
    # For duplicated symbols, append the gene_id to make them unique
    unique_symbols <- symbols
    unique_symbols[duplicated_mask] <- paste0(symbols[duplicated_mask], "_", gene_ids[duplicated_mask])
    
    return(unique_symbols)
  } else {
    return(symbols)
  }
}

# This function would replace the problematic section around line 178-180:
# Instead of:
# rownames(all.reads_c) <- all.reads$Symbol
#
# Use:
# unique_symbols <- make_unique_symbols(all.reads$Symbol, all.reads$gene_ids)
# rownames(all.reads_c) <- unique_symbols

cat("Alternative fix created - this makes gene symbols unique instead of using gene_ids\n")
cat("To use this approach, modify deseq2_norm_qc.R around line 180\n")