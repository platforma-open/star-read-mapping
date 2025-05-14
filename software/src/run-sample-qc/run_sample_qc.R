#!/usr/bin/env Rscript

# This script was tested with R version 4.3.3 and DESeq2 version 1.42.0
# Modified to handle single/no sample inputs without erroring and to produce

# Load libraries
library(DESeq2)
library(reshape2)

# Function to write outputs when processing cannot proceed meaningfully or for empty inputs
writeMinimalOutputs <- function(output_folder, sample_names = character(0)) {
  if (!dir.exists(output_folder)) { dir.create(output_folder, recursive = TRUE) }
  
  # PCA results: Headers "Sample,PC1,PC2"
  if (length(sample_names) > 0) {
    pca_out_df <- data.frame(Sample = sample_names, PC1 = NA_real_, PC2 = NA_real_)
  } else { 
    pca_out_df <- data.frame(Sample = character(0), PC1 = numeric(0), PC2 = numeric(0))
  }
  write.csv(pca_out_df, file = paste0(output_folder, "/pca_results.csv"), row.names = FALSE)

  # Sample distances: Header "Sample" and potentially other sample names if known, but 0 rows.
  # For a truly minimal file with just the "Sample" header:
  dist_minimal_df <- data.frame(Sample=character(0))
  # If sample_names were available and we wanted columns for them (e.g. Sample, S1, S2):
  # if(length(sample_names) > 0) {
  #   for(s_name in sample_names) {
  #     dist_minimal_df[[s_name]] <- numeric(0)
  #   }
  # }
  write.csv(dist_minimal_df, file = paste0(output_folder, "/sample_distances.csv"), row.names = FALSE)
  
  cat("Minimal PCA results and sample distances files created.\n")
}

# Function to process RNA-Seq data
processRawCounts <- function(raw_counts_path, output_folder) {

  # Load counts long format
  countsDataLong <- read.csv(raw_counts_path)
  
  if (nrow(countsDataLong) == 0) {
    cat("Input file is empty. No processing will be done.\n")
    writeMinimalOutputs(output_folder)
    return()
  }

  sample_col <- names(countsDataLong)[1]
  gene_col <- names(countsDataLong)[2]
  count_col <- names(countsDataLong)[3]

  countsDataWide <- dcast(countsDataLong, formula = paste(gene_col, "~", sample_col), value.var = count_col)
  
  if (any(duplicated(countsDataWide[[gene_col]]))) {
    cat("Error: Duplicate gene identifiers found: ", paste(unique(countsDataWide[[gene_col]][duplicated(countsDataWide[[gene_col]])]), collapse=", "), "\n")
    cat("Aborting processing.\n")
    potential_sample_names <- if(ncol(countsDataWide) > 1 && gene_col %in% names(countsDataWide)) colnames(countsDataWide)[!colnames(countsDataWide) %in% gene_col] else character(0)
    if (gene_col %in% names(countsDataWide) && ncol(countsDataWide) == 1 && nrow(countsDataWide) > 0) { 
        potential_sample_names <- character(0)
    }
    writeMinimalOutputs(output_folder, potential_sample_names) 
    return()
  }
  if(all(is.na(countsDataWide[[gene_col]])) || length(countsDataWide[[gene_col]]) == 0) {
    cat("Error: Gene identifier column is empty or all NA. Aborting processing.\n")
    writeMinimalOutputs(output_folder)
    return()
  }
  rownames(countsDataWide) <- countsDataWide[[gene_col]]
  countsDataWide[[gene_col]] <- NULL 
  
  current_rownames <- rownames(countsDataWide)
  if (ncol(countsDataWide) == 0 && nrow(countsDataLong) > 0 && length(unique(countsDataLong[[sample_col]]))==0) {
     cat("No sample columns created after reshaping. Check input format.\n")
     writeMinimalOutputs(output_folder)
     return()
  } else if (ncol(countsDataWide) == 0) {
     cat("Count matrix is empty after initial processing. Aborting.\n")
     writeMinimalOutputs(output_folder)
     return()
  }

  countsDataWide <- as.data.frame(lapply(countsDataWide, function(x) as.numeric(as.character(x))))
  rownames(countsDataWide) <- current_rownames

  num_samples <- ncol(countsDataWide)
  num_genes <- nrow(countsDataWide) 
  current_sample_names <- colnames(countsDataWide)

  if (num_samples == 0) {
    cat("No samples found after reshaping data. Aborting processing.\n")
    writeMinimalOutputs(output_folder) 
    return()
  }

  transformedDataForPcaAndDist <- NULL

  if (num_samples > 1 && num_genes > 0) {
    cat("Multiple samples detected. Proceeding with DESeq2 VST.\n")
    colData <- data.frame(row.names = current_sample_names)
    countsDataWideInteger <- round(countsDataWide) 
    if (is.data.frame(countsDataWideInteger)) {
        countsDataWideInteger <- as.matrix(countsDataWideInteger)
    }
    storage.mode(countsDataWideInteger) <- "integer"
    dds <- DESeqDataSetFromMatrix(countData = countsDataWideInteger, 
                                  colData = DataFrame(condition = rep("none", num_samples)), 
                                  design = ~ 1)
    if (num_genes < 1000) {
      cat("Dataset has less than 1000 genes. Using varianceStabilizingTransformation.\n")
      vstData <- varianceStabilizingTransformation(dds, blind = TRUE)
    } else {
      cat("Dataset has 1000 or more genes. Using vst.\n")
      vstData <- vst(dds, blind = TRUE)
    }
    transformedDataForPcaAndDist <- assay(vstData)
  } else if (num_samples == 1 && num_genes > 0) { 
    cat("Single sample detected. Using log2(counts + 1) transformation.\n")
    transformedDataForPcaAndDist <- log2(countsDataWide + 1)
  } else { 
    cat("Not enough data for reliable transformation. PCA and distances might be trivial or NA.\n")
  }
  
  final_pca_df <- data.frame(Sample = current_sample_names, PC1 = NA_real_, PC2 = NA_real_)

  if (num_genes > 0 && num_samples > 0 && !is.null(transformedDataForPcaAndDist)) {
    scale_pca <- if (num_samples > 1 && num_genes > 1) TRUE else FALSE 
    cat(paste("Performing PCA with center=TRUE and scale=", scale_pca, "\n"))
    data_for_prcomp <- t(transformedDataForPcaAndDist)
    if (!is.matrix(data_for_prcomp)) {
      data_for_prcomp <- as.matrix(data_for_prcomp)
    }
    pcaRes <- NULL
    tryCatch({
        if(scale_pca && any(apply(data_for_prcomp, 2, function(col) var(col, na.rm=TRUE) < 1e-6))) {
            cat("Warning: Some genes have near-zero variance. Disabling scaling for PCA.\n")
            scale_pca <- FALSE
        }
        pcaRes <- prcomp(data_for_prcomp, center = TRUE, scale. = scale_pca)
    }, error = function(e) {
        cat("Error during prcomp: ", e$message, "\nPCA results will use NAs where appropriate.\n")
    })

    if (!is.null(pcaRes)) {
        pca_scores_raw <- as.data.frame(pcaRes$x)
        if (num_samples > 1) {
            final_pca_df <- data.frame(Sample = rownames(pca_scores_raw))
            final_pca_df <- cbind(final_pca_df, pca_scores_raw)
        } else { 
            final_pca_df$Sample[1] <- rownames(pca_scores_raw)[1]
            if ("PC1" %in% colnames(pca_scores_raw)) {
                final_pca_df$PC1[1] <- pca_scores_raw[1, "PC1"]
            }
            if ("PC2" %in% colnames(pca_scores_raw)) { 
                final_pca_df$PC2[1] <- pca_scores_raw[1, "PC2"]
            } else if (num_samples == 1) { 
                 final_pca_df$PC2[1] <- NA_real_
            }
        }
    }
  } else {
      cat("Skipping PCA calculation (no genes, no samples, or transformation failed).\n")
  }
  write.csv(final_pca_df, file = paste0(output_folder, "/pca_results.csv"), row.names = FALSE)
  
  # Sample Distances: "Sample" as first col, then sample names as other col headers.
  dist_df_to_write <- data.frame(Sample=character(0)) # Default for empty/error

  if (num_genes > 0 && !is.null(transformedDataForPcaAndDist)) {
    if (num_samples >= 2) {
      cat("Calculating sample distances for multiple (>=2) samples.\n")
      dist_matrix_temp <- NULL
      tryCatch({
        dist_val <- dist(t(transformedDataForPcaAndDist)) 
        dist_matrix_temp <- as.matrix(dist_val)
      }, error = function(e){
        cat("Error calculating distance matrix: ", e$message, "\nDistances will be minimal.\n")
      })
      if(!is.null(dist_matrix_temp)){
        dist_df_to_write <- as.data.frame(dist_matrix_temp)
        dist_df_to_write <- cbind(Sample = rownames(dist_df_to_write), dist_df_to_write)
      } else { # Error occurred, use current sample names for an empty structure
        dist_df_to_write <- data.frame(Sample=current_sample_names)
        for(s_name in current_sample_names) { dist_df_to_write[[s_name]] <- NA_real_ }
      }
    } else if (num_samples == 1) {
      cat("Processing for a single sample. Distance matrix will be 0 (self-distance).\n")
      sample_name <- current_sample_names[1]
      dist_df_to_write <- data.frame(Sample = sample_name)
      dist_df_to_write[[sample_name]] <- 0 
    }
  } else {
    cat("Skipping distance calculation (no genes, or transformation failed).\n")
    # If current_sample_names exist, create an empty structure with them
    if(length(current_sample_names) > 0) {
        dist_df_to_write <- data.frame(Sample=current_sample_names)
        # Add NA columns for each sample name if we want that structure
        # for(s_name in current_sample_names) { dist_df_to_write[[s_name]] <- NA_real_ }
    } # Else dist_df_to_write remains data.frame(Sample=character(0))
  }
  write.csv(dist_df_to_write, file = paste0(output_folder, "/sample_distances.csv"), row.names = FALSE)

  cat("Processing complete. PCA results and sample distances saved.\n")
}

# Main script body
args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 2) {
  raw_counts_path <- args[1]
  output_folder <- args[2]
  
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    cat("Created output folder:", output_folder, "\n")
  }
  processRawCounts(raw_counts_path, output_folder)
} else {
  cat("Usage: Rscript run_sample_qc.R <raw_counts_path> <output_folder>\n")
}
