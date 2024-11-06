#!/usr/bin/env Rscript

# This script was tested with R version 4.3.3 and DESeq2 version 1.42.0

# Load required libraries
library(DESeq2)

# Function to process RNA-Seq data with VST without detailed sample information
processRawCounts <- function(raw_counts_path, output_folder) {
  # Load raw counts data
  countsData <- read.csv(raw_counts_path, row.names = 1)
  
  # Create a minimal colData with no specific sample information
  colData <- data.frame(row.names = colnames(countsData))
  dds <- DESeqDataSetFromMatrix(countData = countsData, colData = DataFrame(condition = rep("none", ncol(countsData))), design = ~ 1)
  
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
  write.csv(distMatrix, file = paste0(output_folder, "/sampleDistances.csv"))

  distMatrixDf <- as.data.frame(distMatrix)
  distMatrixDf$Sample <- rownames(distMatrixDf)
  distMatrixDf <- distMatrixDf[, c("Sample", names(distMatrixDf)[-ncol(distMatrixDf)])]
  # write.csv(distMatrixDf, file = paste0(output_folder, "/sample_distances.csv"), row.names = FALSE)

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

