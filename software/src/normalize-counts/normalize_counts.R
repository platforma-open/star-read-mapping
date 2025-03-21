#!/usr/bin/env Rscript

# Load libraries
suppressMessages(library(DESeq2))
suppressMessages(library(reshape2))
suppressMessages(library(AnnotationDbi))

# Function to dynamically load species-specific annotation package
load_annotation_package <- function(species) {
  # Map species names to annotation packages
  species_to_package <- list(
    "homo-sapiens" = "org.Hs.eg.db",
    "mus-musculus" = "org.Mm.eg.db",
    "rattus-norvegicus" = "org.Rn.eg.db",
    "danio-rerio" = "org.Dr.eg.db",
    "drosophila-melanogaster" = "org.Dm.eg.db",
    "arabidopsis-thaliana" = "org.At.tair.db",
    "saccharomyces-cerevisiae" = "org.Sc.sgd.db",
    "caenorhabditis-elegans" = "org.Ce.eg.db",
    "gallus-gallus" = "org.Gg.eg.db",
    "bos-taurus" = "org.Bt.eg.db",
    "sus-scrofa" = "org.Ss.eg.db",
    "test-species" = "org.Mm.eg.db"
  )

  # Check if the species name is valid
  if (!(species %in% names(species_to_package))) {
    stop("Unsupported species name. Supported species are: ",
         paste(names(species_to_package), collapse = ", "))
  }

  # Load the appropriate package
  annotation_package <- species_to_package[[species]]
  suppressMessages(library(annotation_package, character.only = TRUE))
  return(annotation_package)
}

# Annotate only matched IDs and keep unmatched IDs in the final output
process_RNASeq <- function(raw_counts_path, output_folder, species) {

  # Load species-specific annotation package
  annotation_package <- load_annotation_package(species)

  # Read raw counts
  countsDataLong <- read.csv(raw_counts_path)

  sample_col <- names(countsDataLong)[1]
  gene_col <- names(countsDataLong)[2]
  count_col <- names(countsDataLong)[3]

  # Reshape to classic count matrix format (sample as col, gene as row)
  countsDataWide <- dcast(countsDataLong,
                          formula = paste(gene_col, "~", sample_col),
                          value.var = count_col)
  rownames(countsDataWide) <- countsDataWide[[gene_col]]
  countsDataWide[[gene_col]] <- NULL

  # Create a minimal colData with no specific sample information
  colData <- data.frame(row.names = colnames(countsDataWide))
  dds <- DESeqDataSetFromMatrix(countData = countsDataWide,
                                colData = DataFrame(condition = rep("none", ncol(countsDataWide))),
                                design = ~ 1)
  dds <- DESeq(dds)

  # Extract normalized counts
  normalized_counts <- counts(dds, normalized = TRUE)

  # Strip version numbers from ENSEMBL keys
  rownames(normalized_counts) <- sub("\\.\\d+$", "",
                                     rownames(normalized_counts))

  # Define Ensembl IDs from the count matrix
  ensembl_ids <- rownames(normalized_counts)  # Assuming normalized_counts has Ensembl IDs as rownames

  # Get all valid Ensembl IDs from the annotation database
  valid_ensembl_ids <- keys(get(annotation_package), keytype = "ENSEMBL")

  # Find matching and non-matching IDs
  matched_ids <- ensembl_ids[ensembl_ids %in% valid_ensembl_ids]
  # unmatched_ids <- ensembl_ids[!ensembl_ids %in% valid_ensembl_ids]

  # Dynamically set the column and key type based on the species
  if (species == "saccharomyces-cerevisiae") {
    column_to_map <- "COMMON"  # For yeast, use COMMON instead of SYMBOL
    key_type <- "ENSEMBL"
  } else if (species == "arabidopsis-thaliana") {
    column_to_map <- "SYMBOL"
    key_type <- "TAIR"         # For Arabidopsis, use TAIR IDs
  } else {
    column_to_map <- "SYMBOL"
    key_type <- "ENSEMBL"      # Default to ENSEMBL for other species
  }

  # Annotate only matched Ensembl IDs
  matched_symbols <- mapIds(
    get(annotation_package),
    keys = matched_ids,
    column = column_to_map,
    keytype = key_type,
    multiVals = "first"
  )

  # Create a SYMBOL column for the entire dataset
  gene_symbols <- sapply(ensembl_ids, function(id) {
    if (id %in% names(matched_symbols)) {
      return(matched_symbols[[id]])
    } else {
      return(NA)  # Leave unmatched IDs with NA
    }
  })

  # Convert row names to a column for melting
  normalized_counts_df <- as.data.frame(normalized_counts)
  normalized_counts_df$Geneid <- rownames(normalized_counts_df)
  normalized_counts_df$SYMBOL <- gene_symbols

  # Convert to long format for output
  normalized_counts_long <- melt(
    normalized_counts_df,
    id.vars = c("Geneid", "SYMBOL"),
    variable.name = "Sample",
    value.name = "NormCounts"
  )

  # Reorder columns
  normalized_counts_long <- normalized_counts_long[, c("Sample", "Geneid",
                                                       "SYMBOL", "NormCounts")]

  # Write normalized counts to CSV
  normalized_counts_path <- file.path(output_folder, "normalized_counts.csv")
  write.csv(normalized_counts_long, file = normalized_counts_path,
            row.names = FALSE)

  # # Log-transform the normalized counts
  # log_normalized_counts <- log2(normalized_counts + 1)

  # # Write log-transformed normalized counts to CSV
  # log_normalized_counts_path <- file.path(output_folder, "log_normalized_counts.csv")
  # write.csv(log_normalized_counts, file = log_normalized_counts_path, row.names = TRUE)

  cat("Normalized counts have been saved to:", normalized_counts_path, "\n")
  #cat("Log transformed normalized counts have been saved to:", log_normalized_counts_path, "\n")
}

# Main execution block
args <- commandArgs(trailingOnly = TRUE)
print(args)

# Check for correct number of arguments
if (length(args) != 3) {
  cat("Usage: Rscript normalize_counts.R <raw_counts_path> <output_folder> <species>\n")
  quit(status = 1)
}

# Assign arguments
raw_counts_path <- args[1]
output_folder <- args[2]
species <- args[3]

# Process RNASeq data
process_RNASeq(raw_counts_path, output_folder, species)
