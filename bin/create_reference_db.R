#!/usr/bin/env Rscript

#
# LMDseq tximport script - Enhanced with comprehensive annotations
# Created by Pierluigi Di Chiaro
# Licensed under MIT License - see LICENSE file for details
#

# Load required libraries
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("tximport")
  library("rtracklayer")
  library("GenomicFeatures")
})

# Helper function to populate gene annotations
populate_gene_annotations <- function(GENES, ii, Anno_Txdb, ww) {
    if("tx_id" %in% names(Anno_Txdb)) GENES$tx_id[ii] <- as.character(Anno_Txdb$tx_id[ww])
    if("seqnames" %in% names(Anno_Txdb)) GENES$seqnames[ii] <- as.character(Anno_Txdb$seqnames[ww])
    if("start" %in% names(Anno_Txdb)) GENES$start[ii] <- as.character(Anno_Txdb$start[ww])
    if("end" %in% names(Anno_Txdb)) GENES$end[ii] <- as.character(Anno_Txdb$end[ww])
    if("width" %in% names(Anno_Txdb)) GENES$width[ii] <- as.character(Anno_Txdb$width[ww])
    if("strand" %in% names(Anno_Txdb)) GENES$strand[ii] <- as.character(Anno_Txdb$strand[ww])
    if("LOC" %in% names(Anno_Txdb)) GENES$LOC[ii] <- as.character(Anno_Txdb$LOC[ww])
    if("type_of_gene" %in% names(Anno_Txdb)) GENES$type_of_gene[ii] <- as.character(Anno_Txdb$type_of_gene[ww])
    if("Symbol" %in% names(Anno_Txdb)) GENES$Symbol[ii] <- as.character(Anno_Txdb$Symbol[ww])
    if("HGNC" %in% names(Anno_Txdb)) GENES$HGNC[ii] <- as.character(Anno_Txdb$HGNC[ww])
    if("Name" %in% names(Anno_Txdb)) GENES$Name[ii] <- as.character(Anno_Txdb$Name[ww])
    if("Entrez.Gene.ID" %in% names(Anno_Txdb)) GENES$Entrez.Gene.ID[ii] <- as.character(Anno_Txdb$Entrez.Gene.ID[ww])
    if("cDNA_LENGTH" %in% names(Anno_Txdb)) GENES$cDNA_LENGTH[ii] <- as.character(Anno_Txdb$cDNA_LENGTH[ww])
    
    return(GENES)
}

# Command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Space-separated list of kallisto output directories", metavar="character"),
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="GTF annotation file", metavar="character"),
  make_option(c("-r", "--reference"), type="character", default=NULL,
              help="Reference annotation file with comprehensive gene annotations", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output file name", metavar="character")
)

optparser <- OptionParser(option_list=option_list)
opt <- parse_args(optparser)

# Validate required arguments
if (is.null(opt$input)) {
  print_help(optparser)
  stop("Input directories must be provided", call.=FALSE)
}
if (is.null(opt$gtf)) {
  print_help(optparser)
  stop("GTF file must be provided", call.=FALSE)
}
if (is.null(opt$output)) {
  print_help(optparser)
  stop("Output file must be provided", call.=FALSE)
}

# Parse input parameters
input_dirs <- trimws(strsplit(opt$input, " ")[[1]])
gtf_file <- opt$gtf
ref_file <- opt$reference
output_file <- opt$output

cat("Input directories:", paste(input_dirs, collapse=", "), "\n")
cat("GTF file:", gtf_file, "\n")
cat("Reference file:", ifelse(is.null(ref_file), "None", ref_file), "\n")
cat("Output file:", output_file, "\n")

# Validate inputs
for(dir in input_dirs) {
    if(!dir.exists(dir)) {
        stop("Directory not found: ", dir, call.=FALSE)
    }
    abundance_file_h5 <- file.path(dir, "abundance.h5")
    abundance_file_tsv <- file.path(dir, "abundance.tsv")
    if(!file.exists(abundance_file_h5) && !file.exists(abundance_file_tsv)) {
        stop("abundance.h5 or abundance.tsv not found in: ", dir, call.=FALSE)
    }
}

if(!file.exists(gtf_file)) {
    stop("GTF file not found: ", gtf_file, call.=FALSE)
}

## ---- Associate transcripts to gene IDs by tximport ---------
cat("Associating transcripts to gene IDs\n")

# Read GTF and create transcript-to-gene mapping
gtf_data <- import(gtf_file)
gtf_df <- as.data.frame(gtf_data)

# Extract transcript and gene information
gencode_transcript <- dplyr::filter(gtf_df, type=="transcript")
gencode_gene <- dplyr::filter(gtf_df, type=="gene")

if(nrow(gencode_transcript) == 0) {
    stop("No transcripts found in GTF file", call.=FALSE)
}

tx2gene <- dplyr::select(gencode_transcript, transcript_id, gene_id)
tx2gene <- tx2gene[complete.cases(tx2gene), ]

cat("Found", nrow(tx2gene), "transcript-to-gene mappings\n")

# Prepare file paths for tximport
samples <- basename(input_dirs)

# Check file type and prepare accordingly
use_h5 <- all(file.exists(file.path(input_dirs, "abundance.h5")))
use_tsv <- all(file.exists(file.path(input_dirs, "abundance.tsv")))

if(use_h5) {
    abundance_files <- file.path(input_dirs, "abundance.h5")
    file_type <- "kallisto"
    cat("Using H5 files for tximport...\n")
} else if(use_tsv) {
    abundance_files <- file.path(input_dirs, "abundance.tsv")
    file_type <- "kallisto"
    cat("Using TSV files for tximport...\n")
} else {
    stop("No consistent abundance files found across all samples", call.=FALSE)
}

names(abundance_files) <- samples

# Run tximport
cat("Running tximport...\n")
txi.kallisto <- tximport(abundance_files, type = file_type, tx2gene = tx2gene)
all.reads <- as.data.frame(txi.kallisto$counts)
names(all.reads) <- samples

cat("Imported expression data for", nrow(all.reads), "genes and", ncol(all.reads), "samples\n")

## ---- Create the reference on which we will record the comprehensive annotations ---------
cat("Creating comprehensive transcript database\n")

# Initialize comprehensive gene annotation dataframe
GENES <- as.data.frame(matrix(NA, nrow=nrow(all.reads)))
rownames(GENES) <- rownames(all.reads)
colnames(GENES) <- "gene_id"

GENES$gene_id <- as.character(rownames(GENES))
GENES$gene_ids <- gsub(GENES$gene_id, pattern="\\.[0-9]+$|\\.[0-9]+_PAR_Y", replacement="")
GENES$IDs <- NA
GENES$tx_id <- NA
GENES$seqnames <- NA
GENES$start <- NA
GENES$end <- NA
GENES$width <- NA
GENES$strand <- NA
GENES$LOC <- NA
GENES$type_of_gene <- NA
GENES$Symbol <- NA
GENES$HGNC <- NA
GENES$Name <- NA
GENES$Entrez.Gene.ID <- NA
GENES$cDNA_LENGTH <- NA

# Extract gene names from GTF
cat("Extracting gene names from GTF...\n")
for(ii in seq_along(GENES$gene_id)){
  if(ii %% 1000 == 0) cat(ii," of ",length(GENES$gene_id),"\n")
  ww <- which(gencode_gene$gene_id == GENES$gene_id[ii])
  if(length(ww)!=0){
    GENES$IDs[ii] <- as.character(gencode_gene$gene_name[ww])
  }
}

# Handle duplicated gene names
cat("Handling duplicated gene names...\n")
dup <- GENES$IDs[duplicated(GENES$gene_ids)]

for(ii in dup){
  if(!is.na(ii)) {
    w <- which(GENES$IDs==ii)
    if(length(w) == 1){
      cat(ii," not duplicated \n")
    }else{
      cat(ii," is duplicated \n")
      for(jj in seq_along(w)) {
        GENES[w[jj],"IDs"] <- paste0(ii,":",jj)
      }
    }
  }
}

## ---- Import reference file and add comprehensive annotations ---------
if(!is.null(ref_file) && file.exists(ref_file)) {
    cat("Loading comprehensive reference annotations from:", ref_file, "\n")
    
    tryCatch({
        # Read reference annotation file
        Anno_Txdb <- read.delim(ref_file, stringsAsFactors = FALSE)
        
        # Set rownames if there's an 'X' column (common in reference files)
        if("X" %in% names(Anno_Txdb)) {
            rownames(Anno_Txdb) <- Anno_Txdb$X
        }
        
        cat("Reference file loaded with", nrow(Anno_Txdb), "annotations\n")
        cat("Available columns:", paste(names(Anno_Txdb), collapse=", "), "\n")
        
        # Add comprehensive annotations using multiple matching strategies
        cat("Matching genes to reference annotations...\n")
        matched_count <- 0
        
        for(ii in seq_along(GENES$gene_ids)){
            if(ii %% 1000 == 0) cat(ii," of ",length(GENES$gene_ids),"\n")
            
            # Strategy 1: Match by ENSEMBL ID (removing version)
            ww <- NULL
            if("ENSEMBL" %in% names(Anno_Txdb)) {
                ww <- grep(GENES$gene_ids[ii], gsub("\\;ENSG*", "", Anno_Txdb$ENSEMBL))
            }
            
            if(length(ww)==1){
                matched_count <- matched_count + 1
                GENES <- populate_gene_annotations(GENES, ii, Anno_Txdb, ww)
            } else {
                # Strategy 2: Match by ENSEMBL_CHECK if available
                if("ENSEMBL_CHECK" %in% names(Anno_Txdb)) {
                    ww <- grep(GENES$gene_ids[ii], gsub("\\;ENSG*", "", Anno_Txdb$ENSEMBL_CHECK))
                    if(length(ww)==1){
                        matched_count <- matched_count + 1
                        GENES <- populate_gene_annotations(GENES, ii, Anno_Txdb, ww)
                    } else if(length(ww)==0) {
                        # Strategy 3: Match by gene symbol/name
                        if(!is.na(GENES$IDs[ii])) {
                            ww <- which(rownames(Anno_Txdb) == GENES$IDs[ii])
                            if(length(ww)==1){
                                matched_count <- matched_count + 1
                                GENES <- populate_gene_annotations(GENES, ii, Anno_Txdb, ww)
                            }
                        }
                    }
                }
            }
        }
        
        cat("Successfully matched", matched_count, "out of", nrow(GENES), "genes (", 
            round(matched_count/nrow(GENES)*100, 1), "%)\n")
        
    }, error = function(e) {
        cat("Warning: Could not process reference file:", e$message, "\n")
        cat("Proceeding with basic GTF annotations only\n")
    })
} else {
    cat("No reference file provided - using GTF annotations only\n")
}

## ---- Handle duplicate gene symbols ---------
cat("Handling duplicate gene symbols...\n")

# Handle duplicates in Symbol column using the same approach as IDs column
dup_symbols <- GENES$Symbol[duplicated(GENES$Symbol)]
unique_dup_symbols <- unique(dup_symbols[!is.na(dup_symbols)])

if(length(unique_dup_symbols) > 0) {
    cat("Found", length(unique_dup_symbols), "duplicated symbols to fix\n")
    
    for(symbol in unique_dup_symbols){
        w <- which(GENES$Symbol == symbol & !is.na(GENES$Symbol))
        if(length(w) == 1){
            cat(symbol, " not duplicated \n")
        } else {
            cat(symbol, " is duplicated (", length(w), " times)\n")
            # Apply same numbering scheme as IDs: SYMBOL:1, SYMBOL:2, etc.
            for(j in seq_along(w)) {
                GENES$Symbol[w[j]] <- paste0(symbol, ":", j)
            }
        }
    }
} else {
    cat("No duplicate symbols found\n")
}

## ---- Handle duplicate gene symbols using the same approach as IDs ---------
cat("Handling duplicate gene symbols...\n")

# Use the exact same logic as for IDs column
dup <- GENES$Symbol[duplicated(GENES$Symbol) & !is.na(GENES$Symbol)]

for(ii in dup){ 
  w <- which(GENES$Symbol==ii) 
  if(length(w) == 1){ 
    cat(ii," not duplicated \n") 
  }else{ 
    cat(ii," is duplicated \n") 
    GENES[w[1],"Symbol"] <- paste0(ii,":1") 
    GENES[w[2],"Symbol"] <- paste0(ii,":2")
    # Handle cases with more than 2 duplicates
    if(length(w) > 2) {
      for(j in 3:length(w)) {
        GENES[w[j],"Symbol"] <- paste0(ii,":",j)
      }
    }
  } 
}

## ---- Filter and clean the data ---------
cat("Filtering and cleaning data...\n")

# Remove genes with no Symbol annotation
initial_count <- nrow(GENES)
set <- which(is.na(GENES$Symbol))
if(length(set)!=0){
    cat("Removing", length(set), "genes with no Symbol annotation\n")
    GENES <- GENES[-set,]
}

# Filter out genes on the chrY PAR regions
set <- grep("_PAR_Y", GENES$gene_id)
if(length(set)!=0){
    GENES <- GENES[-set,]
    cat("Removed", length(set), "PAR_Y genes\n")
}

# Maintain only genes which are present in both NCBI and ENSEMBL
if(!all(is.na(GENES$Symbol)) && !all(is.na(GENES$IDs))) {
    set <- which(GENES$Symbol == gsub("\\:[0-9]+$", "", GENES$IDs))
    if(length(set)!=0){
        cat("Keeping", length(set), "genes with matching Symbol and gene names (present in both NCBI and ENSEMBL)\n")
        GENES <- GENES[set,]
    } else {
        cat("Warning: No genes found with matching Symbol and IDs - this may indicate annotation issues\n")
    }
}

cat("Final gene count:", nrow(GENES), "genes retained from initial", initial_count, "\n")

## ---- Filter expression data and create final output ---------
cat("Creating final output matrix...\n")

# Filter expression data to match cleaned gene annotations
reads <- all.reads[rownames(GENES), samples, drop=FALSE]

# Convert estimated counts to integers (required for DESeq2 downstream)
cat("Rounding estimated counts to integers...\n")
reads[, samples] <- round(reads[, samples])

RES_DAT_ALL <- cbind(GENES, reads)

# Order by transcript ID if available
if(!all(is.na(RES_DAT_ALL$tx_id))) {
    numeric_tx_ids <- suppressWarnings(as.numeric(RES_DAT_ALL$tx_id))
    if(!all(is.na(numeric_tx_ids))) {
        RES_DAT_ALL <- RES_DAT_ALL[order(numeric_tx_ids, na.last=TRUE),]
    }
}

# Remove tx_id column and reset row names
RES_DAT_ALL$tx_id <- NULL
rownames(RES_DAT_ALL) <- NULL

# Write output
cat("Writing output matrix...\n")
write.table(RES_DAT_ALL, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== COMPREHENSIVE ANNOTATION SUMMARY ===\n")
cat("Final genes:", nrow(RES_DAT_ALL), "\n")
cat("Samples:", length(samples), "\n")
cat("Annotation columns:", ncol(GENES), "\n")

# Summary of annotation completeness
annotation_cols <- c("Symbol", "HGNC", "Name", "Entrez.Gene.ID", "LOC", "type_of_gene")
for(col in annotation_cols) {
    if(col %in% names(RES_DAT_ALL)) {
        complete <- sum(!is.na(RES_DAT_ALL[[col]]))
        pct <- round(complete/nrow(RES_DAT_ALL)*100, 1)
        cat(col, ":", complete, "complete (", pct, "%)\n")
    }
}

# Helper function to populate gene annotations
populate_gene_annotations <- function(GENES, ii, Anno_Txdb, ww) {
    if("tx_id" %in% names(Anno_Txdb)) GENES$tx_id[ii] <- as.character(Anno_Txdb$tx_id[ww])
    if("seqnames" %in% names(Anno_Txdb)) GENES$seqnames[ii] <- as.character(Anno_Txdb$seqnames[ww])
    if("start" %in% names(Anno_Txdb)) GENES$start[ii] <- as.character(Anno_Txdb$start[ww])
    if("end" %in% names(Anno_Txdb)) GENES$end[ii] <- as.character(Anno_Txdb$end[ww])
    if("width" %in% names(Anno_Txdb)) GENES$width[ii] <- as.character(Anno_Txdb$width[ww])
    if("strand" %in% names(Anno_Txdb)) GENES$strand[ii] <- as.character(Anno_Txdb$strand[ww])
    if("LOC" %in% names(Anno_Txdb)) GENES$LOC[ii] <- as.character(Anno_Txdb$LOC[ww])
    if("type_of_gene" %in% names(Anno_Txdb)) GENES$type_of_gene[ii] <- as.character(Anno_Txdb$type_of_gene[ww])
    if("Symbol" %in% names(Anno_Txdb)) GENES$Symbol[ii] <- as.character(Anno_Txdb$Symbol[ww])
    if("HGNC" %in% names(Anno_Txdb)) GENES$HGNC[ii] <- as.character(Anno_Txdb$HGNC[ww])
    if("Name" %in% names(Anno_Txdb)) GENES$Name[ii] <- as.character(Anno_Txdb$Name[ww])
    if("Entrez.Gene.ID" %in% names(Anno_Txdb)) GENES$Entrez.Gene.ID[ii] <- as.character(Anno_Txdb$Entrez.Gene.ID[ww])
    if("cDNA_LENGTH" %in% names(Anno_Txdb)) GENES$cDNA_LENGTH[ii] <- as.character(Anno_Txdb$cDNA_LENGTH[ww])
    
    return(GENES)
}

cat("Output file:", output_file, "\n")
cat("=== tximport completed successfully ===\n")

