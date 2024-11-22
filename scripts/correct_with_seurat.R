#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(Seurat)
  library(dplyr)
})

# Define command-line options
option_list <- list(
  make_option(c("--input_data"), type = "character", help = "Path to input data file"),
  make_option(c("--batch_key"), type = "character", help = "Batch column name"),
  make_option(c("--method"), type = "character", help = "Seurat integration method"),
  make_option(c("--output_path"), type = "character", help = "Path to output file")
)

# Parse command-line arguments
parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

# Validate arguments
if (is.null(opt$input_data) || is.null(opt$batch_key) || is.null(opt$method) || is.null(opt$output_path)) {
  print_help(parser)
  stop("All four arguments (--input_data, --batch_key, --method, --output_path) must be supplied.", call. = FALSE)
}

# Assign arguments to variables
input_file <- opt$input_data
batch_col <- opt$batch_key
seurat_method <- opt$method
output_file <- opt$output_path

# Set future globals size (optional, adjust as needed)
options(future.globals.maxSize = 12000 * 1024^2)

# Read the input parquet data
parquet_data <- as.data.frame(read_parquet(input_file))

# Separate metadata and features
col_names <- names(parquet_data)
metadata_cols <- col_names[grepl("^Metadata_", col_names)]
features_cols <- col_names[!grepl("^Metadata_", col_names)]

features <- parquet_data[, features_cols]
metadata <- parquet_data[, metadata_cols]
batch_info <- parquet_data[[batch_col]]  # Use [[ ]] for column selection

# Identify unique batches
batch_names <- unique(batch_info)

# Split features and metadata by batch
batches <- split(features, batch_info)
meta_batches <- split(metadata, batch_info)

# Initialize list to hold Seurat objects
seurat_lists <- list()

# Iterate over each batch to create Seurat objects
for (i in seq_along(batches)) {
  raw <- batches[[i]] %>% as.data.frame()
  meta <- meta_batches[[i]] %>% as.data.frame()
  
  # Clean column names
  names(raw) <- gsub("_", "-", names(raw))
  
  # Create Seurat object
  obj <- CreateSeuratObject(counts = t(raw), meta.data = meta)
  obj <- SetAssayData(object = obj, layer = "data", new.data = t(raw))
  obj <- SetAssayData(object = obj, layer = "scale.data", new.data = t(raw))
  
  # Run PCA (specify features to skip automatic feature selection)
  obj <- RunPCA(obj, features = colnames(raw), verbose = FALSE)
  
  # Append to the list
  seurat_lists[[i]] <- obj
}

# Integrate data using Seurat
anchor_set <- FindIntegrationAnchors(
  object.list = seurat_lists, 
  reduction = seurat_method, 
  anchor.features = colnames(raw), 
  scale = FALSE,
  verbose = FALSE
)

integrated_obj <- IntegrateData(anchor_set, verbose = FALSE)

# Extract integrated data
int_data <- GetAssayData(object = integrated_obj, assay = "integrated", slot = "data")
corrected <- as.data.frame(t(as.matrix(int_data)))
colnames(corrected) <- gsub("-", "_", colnames(corrected))

# Extract and clean metadata (excluding first three columns)
meta_corrected <- integrated_obj@meta.data[, -(1:3)]

# Combine metadata and corrected data
corrected_df <- cbind(meta_corrected, corrected)

# Write the corrected data to the output parquet file
write_parquet(corrected_df, output_file)
