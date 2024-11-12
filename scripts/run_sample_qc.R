#!/usr/bin/env Rscript

# This script was tested with R version 4.3.3 and DESeq2 version 1.42.0

# Install DESeq2 if necessary
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  cat("Installing DESeq2 from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")
  }
  BiocManager::install("DESeq2")
}

# Install reshape2 if necessary
if (!requireNamespace("reshape2", quietly = TRUE)) {
  cat("reshape2 package not found. Installing reshape2...\n")
  install.packages("reshape2", repos = "https://cloud.r-project.org/")
}

# Load libraries
library(DESeq2)
library(reshape2)

# Function to process RNA-Seq data with VST without detailed sample information
processRawCounts <- function(raw_counts_path, output_folder) {

  # Load counts long format
  countsDataLong <- read.csv(raw_counts_path)
  
  sample_col <- names(countsDataLong)[1]
  gene_col <- names(countsDataLong)[2]
  count_col <- names(countsDataLong)[3]

  # Reshape to classic count matrix format (sample as col, gene as row)
  countsDataWide <- dcast(countsDataLong, formula = paste(gene_col, "~", sample_col), value.var = count_col)
  rownames(countsDataWide) <- countsDataWide[[gene_col]]
  countsDataWide[[gene_col]] <- NULL
  
  # Create a minimal colData with no specific sample information
  colData <- data.frame(row.names = colnames(countsDataWide))
  dds <- DESeqDataSetFromMatrix(countData = countsDataWide, colData = DataFrame(condition = rep("none", ncol(countsDataWide))), design = ~ 1)
  
  # Perform variance stabilizing transformation without fitting to a design
  vstData <- vst(dds, blind = TRUE)
  
  # PCA on the VST data
  pcaRes <- prcomp(t(assay(vstData)))
  pcaData <- as.data.frame(pcaRes$x)
  pcaData$Sample <- rownames(pcaData) # Add Sample column from rownames
  pcaData <- pcaData[, c("Sample", names(pcaData)[-ncol(pcaData)])]  # Reorder to have Sample as first column
  write.csv(pcaData, file = paste0(output_folder, "/pca_results.csv"), row.names = FALSE)
  
  # Calculate sample distances based on VST data
  distMatrix <- dist(t(assay(vstData)))
  distMatrix <- as.matrix(distMatrix)
  distMatrixDf <- as.data.frame(distMatrix)
  distMatrixDf$Sample <- rownames(distMatrixDf)
  distMatrixDf <- distMatrixDf[, c("Sample", names(distMatrixDf)[-ncol(distMatrixDf)])]
  write.csv(distMatrixDf, file = paste0(output_folder, "/sample_distances.csv"), row.names = FALSE)

  cat("PCA results and sample distances have been saved to the specified output folder.\n")
}

# Main script body
args <- commandArgs(trailingOnly = TRUE)

# Expecting two arguments: raw_counts_path and output_folder
if(length(args) == 2) {
  raw_counts_path <- args[1]
  output_folder <- args[2]
  processRawCounts(raw_counts_path, output_folder)
} else {
  cat("Usage: Rscript run_sample_qc.R <raw_counts_path> <output_folder>\n")
}

