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
    "test-species" = "org.Mm.eg.db" # Example: test species maps to mouse
  )

  # Check if the species name is valid
  if (!(species %in% names(species_to_package))) {
    stop("Unsupported species name. Supported species are: ",
         paste(names(species_to_package), collapse = ", "))
  }

  # Load the appropriate package
  annotation_package_name <- species_to_package[[species]]
  cat("Loading annotation package:", annotation_package_name, "\n")
  # Ensure the package is installed, if not, stop with a message
  if (!requireNamespace(annotation_package_name, quietly = TRUE)) {
    stop(paste("Annotation package", annotation_package_name, "is not installed. Please install it first, e.g., BiocManager::install('", annotation_package_name, "')", sep=""))
  }
  suppressMessages(library(annotation_package_name, character.only = TRUE))
  return(annotation_package_name)
}

# Function to write an empty normalized counts file with headers
writeMinimalNormalizedOutput <- function(output_folder) {
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  normalized_counts_path <- file.path(output_folder, "normalized_counts.csv")
  # Define expected headers for an empty file
  headers <- "Sample,Geneid,SYMBOL,NormCounts"
  writeLines(headers, con = normalized_counts_path)
  cat("Minimal normalized_counts.csv file with headers created at:", normalized_counts_path, "\n")
}


# Annotate only matched IDs and keep unmatched IDs in the final output
process_RNASeq <- function(raw_counts_path, output_folder, species) {

  # Create output folder if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    cat("Created output folder:", output_folder, "\n")
  }

  # Load species-specific annotation package
  annotation_package <- NULL # Stores the package *name*
  annotation_db_obj <- NULL    # Stores the loaded package *object* for get()
  tryCatch({
    annotation_package <- load_annotation_package(species) # Gets the package name
    annotation_db_obj <- get(annotation_package)          # Loads the actual db object
  }, error = function(e) {
    cat("Error loading annotation package:", e$message, "\n")
    writeMinimalNormalizedOutput(output_folder)
    stop(e) # Re-throw the error to halt script execution
  })


  # Read raw counts
  countsDataLong <- read.csv(raw_counts_path)
  if (nrow(countsDataLong) == 0) {
    cat("Input file", raw_counts_path, "is empty. No processing will be done.\n")
    writeMinimalNormalizedOutput(output_folder)
    return()
  }

  sample_col <- names(countsDataLong)[1]
  gene_col <- names(countsDataLong)[2]
  count_col <- names(countsDataLong)[3]

  # Reshape to classic count matrix format (sample as col, gene as row)
  countsDataWide <- dcast(countsDataLong,
                          formula = paste(gene_col, "~", sample_col),
                          value.var = count_col)
  
  # Check for issues after dcast
  if (nrow(countsDataWide) == 0 || (ncol(countsDataWide) <= 1) ) { 
      cat("Error: Data reshaping resulted in no gene rows or no sample columns.\n")
      cat("Please check the input file format and column names.\n")
      writeMinimalNormalizedOutput(output_folder)
      return()
  }
  if (any(duplicated(countsDataWide[[gene_col]]))) {
    cat("Error: Duplicate gene identifiers found after reshaping:", paste(unique(countsDataWide[[gene_col]][duplicated(countsDataWide[[gene_col]])]), collapse=", "), "\n")
    cat("Aborting processing.\n")
    writeMinimalNormalizedOutput(output_folder)
    return()
  }
  if(all(is.na(countsDataWide[[gene_col]])) || length(countsDataWide[[gene_col]]) == 0) {
    cat("Error: Gene identifier column is empty or all NA after reshaping. Aborting processing.\n")
    writeMinimalNormalizedOutput(output_folder)
    return()
  }

  rownames(countsDataWide) <- countsDataWide[[gene_col]]
  countsDataWide[[gene_col]] <- NULL # Remove gene ID column

  # Ensure counts are numeric after dcast
  current_rnames <- rownames(countsDataWide)
  current_cnames <- colnames(countsDataWide)
  countsDataWide <- as.data.frame(lapply(countsDataWide, function(x) as.numeric(as.character(x))))
  rownames(countsDataWide) <- current_rnames
  colnames(countsDataWide) <- current_cnames


  num_samples <- ncol(countsDataWide)
  num_genes <- nrow(countsDataWide)

  if (num_samples == 0) {
    cat("No sample columns found in the count matrix. Aborting.\n")
    writeMinimalNormalizedOutput(output_folder)
    return()
  }
  if (num_genes == 0) {
    cat("No gene rows found in the count matrix. Aborting.\n")
    writeMinimalNormalizedOutput(output_folder)
    return()
  }
  
  # Ensure countData is integer matrix for DESeqDataSetFromMatrix (though not strictly needed for CPM)
  countsDataWideInteger <- round(countsDataWide)
  if(is.data.frame(countsDataWideInteger)) {
    countsDataWideInteger <- as.matrix(countsDataWideInteger)
  }
  storage.mode(countsDataWideInteger) <- "integer"

  
  normalized_counts <- NULL # Initialize
  if (num_samples > 1) {
    cat("Multiple samples detected. Running DESeq() for normalization.\n")
    # Create a minimal colData with no specific sample information for DESeq2
    colData <- data.frame(row.names = colnames(countsDataWideInteger))
    dds <- DESeqDataSetFromMatrix(countData = countsDataWideInteger,
                                  colData = DataFrame(condition = rep("none", num_samples)),
                                  design = ~ 1)
    dds <- DESeq(dds) 
    normalized_counts <- counts(dds, normalized = TRUE)
  } else { # Single sample case: Use CPM
    cat("Single sample detected. Calculating Counts Per Million (CPM).\n")
    # countsDataWideInteger is a matrix with one column for the single sample
    # Ensure it's treated as numeric for colSums
    single_sample_counts <- as.numeric(countsDataWideInteger[,1])
    total_counts <- sum(single_sample_counts)
    if (total_counts > 0) {
      cpm_values <- (single_sample_counts / total_counts) * 1e6
    } else {
      cat("Warning: Total counts for the single sample is 0. CPM values will be 0.\n")
      cpm_values <- rep(0, length(single_sample_counts))
    }
    # Reshape CPM values back into a matrix/data.frame format consistent with multi-sample output
    normalized_counts <- as.data.frame(matrix(cpm_values, ncol = 1))
    colnames(normalized_counts) <- current_cnames[1] # Assign the original sample name
    rownames(normalized_counts) <- current_rnames   # Assign original gene names
  }

  if (is.null(normalized_counts)) {
      cat("Error: Normalized counts could not be generated. Aborting.\n")
      writeMinimalNormalizedOutput(output_folder)
      return()
  }

  # Strip version numbers from ENSEMBL keys in rownames
  # Ensure normalized_counts is a matrix or data.frame before using rownames
  if (!is.matrix(normalized_counts) && !is.data.frame(normalized_counts)) {
      cat("Warning: normalized_counts is not a matrix or data frame. Skipping version stripping.\n")
      ensembl_ids <- character(0) 
  } else {
      rownames_norm_counts <- rownames(normalized_counts)
      if(is.null(rownames_norm_counts)){
          cat("Warning: Rownames of normalized_counts are NULL. Cannot strip versions or get Ensembl IDs.\n")
          ensembl_ids <- character(0)
      } else {
          # Ensure normalized_counts is a matrix for rownames<- assignment if it's a single-column data frame
          if(is.data.frame(normalized_counts) && ncol(normalized_counts) == 1) {
              temp_matrix <- as.matrix(normalized_counts) # Convert to matrix
              rownames(temp_matrix) <- sub("\\.\\d+$", "", rownames_norm_counts)
              normalized_counts <- as.data.frame(temp_matrix) # Convert back to data.frame
          } else { # If already matrix or multi-column data.frame
              rownames(normalized_counts) <- sub("\\.\\d+$", "", rownames_norm_counts)
          }
          ensembl_ids <- rownames(normalized_counts)
      }
  }

  # --- Annotation part reverted to user's original logic ---
  gene_symbols <- rep(NA_character_, length(ensembl_ids))
  names(gene_symbols) <- ensembl_ids

  if (length(ensembl_ids) > 0 && !is.null(annotation_db_obj)) {
      cat("Starting gene ID annotation...\n")
      
      valid_db_ids <- character(0)
      tryCatch({
          valid_db_ids <- keys(annotation_db_obj, keytype = "ENSEMBL")
          cat(length(valid_db_ids), "valid ENSEMBL keys found in", annotation_package, "\n")
      }, error = function(e) {
          cat("Warning: Could not retrieve ENSEMBL keys from", annotation_package, ". Error:", e$message, "\n")
          cat("Proceeding without pre-filtering against database ENSEMBL keys.\n")
      })

      matched_ids <- ensembl_ids[ensembl_ids %in% valid_db_ids]
      if(length(valid_db_ids) == 0 && length(ensembl_ids) > 0) {
          cat("No valid ENSEMBL keys retrieved from DB or an error occurred. Attempting to map all input IDs.\n")
          matched_ids <- ensembl_ids 
      }
      cat(length(matched_ids), "IDs from data matched with ENSEMBL keys in the database.\n")

      column_to_map <- "SYMBOL"
      key_type_for_mapids <- "ENSEMBL" 

      if (species == "saccharomyces-cerevisiae") {
        column_to_map <- "COMMON"
      } else if (species == "arabidopsis-thaliana") {
        column_to_map <- "SYMBOL"
        key_type_for_mapids <- "TAIR" 
      }

      cat("Using keytype '", key_type_for_mapids, "' for mapIds to retrieve column '", column_to_map, "'.\n")

      if (length(matched_ids) > 0) {
          mapped_values <- tryCatch({
              mapIds(
                annotation_db_obj,
                keys = matched_ids, 
                column = column_to_map,
                keytype = key_type_for_mapids, 
                multiVals = "first"
              )
          }, error = function(e){
              cat("Warning: mapIds failed with error: ", e$message, "\nGene symbols will be NA for these IDs.\n")
              return(NULL)
          })
          
          if(!is.null(mapped_values)){
              for(id in ensembl_ids) { 
                  if (id %in% names(mapped_values)) { 
                      gene_symbols[id] <- mapped_values[[id]]
                  } 
              }
          }
      } else {
          cat("No matched IDs to annotate.\n")
      }
  } else if (is.null(annotation_db_obj)) {
      cat("Annotation database object was not loaded. Skipping gene symbol annotation.\n")
  } else {
      cat("No Ensembl IDs found to annotate.\n")
  }
  # --- End of reverted annotation logic ---


  # Ensure normalized_counts is a data.frame for consistent column binding
  if(is.matrix(normalized_counts)) {
    normalized_counts_df <- as.data.frame(normalized_counts)
  } else if (is.data.frame(normalized_counts)) { # handles multi-sample and single-sample (already df)
    normalized_counts_df <- normalized_counts
  } else { # Should not happen if logic above is correct
    cat("Error: normalized_counts is in an unexpected format. Attempting to convert.\n")
    normalized_counts_df <- data.frame(normalized_counts)
    # Try to assign colnames/rownames if possible, though this is a fallback
    if(ncol(normalized_counts_df) == length(current_cnames)) colnames(normalized_counts_df) <- current_cnames
    if(nrow(normalized_counts_df) == length(current_rnames)) rownames(normalized_counts_df) <- current_rnames
  }


  normalized_counts_df$Geneid <- rownames(normalized_counts_df) 
  normalized_counts_df$SYMBOL <- gene_symbols[normalized_counts_df$Geneid] 

  id_variables = c("Geneid", "SYMBOL")
  if (!all(id_variables %in% names(normalized_counts_df))) {
      if(!"SYMBOL" %in% names(normalized_counts_df)) normalized_counts_df$SYMBOL <- NA_character_
      if(!"Geneid" %in% names(normalized_counts_df) && !is.null(rownames(normalized_counts_df))) {
          normalized_counts_df$Geneid <- rownames(normalized_counts_df)
      } else if (!"Geneid" %in% names(normalized_counts_df)) {
          stop("Critical error: Geneid column is missing and cannot be recovered from rownames before melt operation.")
      }
  }
  
  measure_vars = colnames(normalized_counts_df)[!colnames(normalized_counts_df) %in% id_variables]
  if(length(measure_vars) == 0 && num_samples > 0){ 
      cat("Warning: No sample data columns found in normalized_counts_df for melting. This might indicate an issue.\n")
      if(length(current_cnames) > 0 && !any(current_cnames %in% id_variables)) {
          measure_vars <- current_cnames
          cat("Using original column names as measure_vars:", paste(measure_vars, collapse=", "), "\n")
      }
  }

  if(length(measure_vars) == 0){
      cat("No sample data columns to melt. Output CSV will be minimal or reflect only Geneid/SYMBOL.\n")
      if (nrow(normalized_counts_df) > 0) {
        final_output_df <- data.frame(
            Sample = NA_character_, # Or perhaps the single sample name if num_samples == 1
            Geneid = normalized_counts_df$Geneid,
            SYMBOL = normalized_counts_df$SYMBOL,
            NormCounts = NA_real_ # Or the actual counts if they exist for the single sample
        )
        if(num_samples == 1 && !is.null(normalized_counts_df[[current_cnames[1]]])) {
            final_output_df$Sample <- current_cnames[1]
            final_output_df$NormCounts <- normalized_counts_df[[current_cnames[1]]]
        }

      } else {
        final_output_df <- data.frame(Sample=character(), Geneid=character(), SYMBOL=character(), NormCounts=numeric())
      }
  } else {
      normalized_counts_long <- melt(
        normalized_counts_df,
        id.vars = id_variables,
        measure.vars = measure_vars, 
        variable.name = "Sample",
        value.name = "NormCounts"
      )
      final_output_df <- normalized_counts_long[, c("Sample", "Geneid", "SYMBOL", "NormCounts")]
  }

  normalized_counts_path <- file.path(output_folder, "normalized_counts.csv")
  write.csv(final_output_df, file = normalized_counts_path, row.names = FALSE, na = "") 

  cat("Normalized counts have been saved to:", normalized_counts_path, "\n")
}

# Main execution block
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  cat("Usage: Rscript normalize_counts.R <raw_counts_path> <output_folder> <species>\n")
  quit(status = 1)
}

raw_counts_path <- args[1]
output_folder <- args[2]
species <- args[3]

process_RNASeq(raw_counts_path, output_folder, species)
