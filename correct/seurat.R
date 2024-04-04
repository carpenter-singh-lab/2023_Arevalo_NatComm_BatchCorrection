args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  cat("Usage: Rscript seurat.R input_file output_file batch_col\n")
  quit()
}
input_file <- args[1]
output_file <- args[2]
batch_col <- args[3]
seurat_method <- args[4]

library(arrow)
library(Seurat)
library(dplyr)

options(future.globals.maxSize = 10000 * 1024^2)
parquet_data <- as.data.frame(read_parquet(input_file))
col_names <- names(parquet_data)
metadata_cols <- col_names[grepl("^Metadata_", col_names)]
features_cols <- col_names[!grepl("^Metadata_", col_names)]
features <- parquet_data[, features_cols]
metadata <- parquet_data[, metadata_cols]
batch_info <- parquet_data[, batch_col]

batch_names <- unique(batch_info)

batches <- split(features, batch_info)
meta_batches <- split(metadata, batch_info)

seurat_lists <- list()
for (i in c(1:length(batches))) {
  raw <- batches[[i]] %>% as.data.frame()
  meta <- meta_batches[[i]] %>% as.data.frame()
  
  # Create a Seurat object from the batch
  names(raw) <- gsub("_","-",names(raw))
  obj <- CreateSeuratObject(counts = t(raw), meta.data = meta)
  obj[["RNA"]]$data <- t(raw)
  obj[["RNA"]]$scale.data <- t(raw)
  
  # Specifying features enables skipping automatic trigger of "FindVariableFeatures" function
  obj <- RunPCA(obj, features = colnames(raw), verbose = FALSE)
  
  seurat_lists[i] <- obj
}

# Integrate using Seurat
anchor_set <- FindIntegrationAnchors(object.list = seurat_lists, 
                                     reduction = seurat_method, 
                                     anchor.features = colnames(raw), 
                                     scale = FALSE,
                                     verbose = FALSE)
integrated_obj = IntegrateData(anchor_set, verbose = FALSE)

# Extract metadata and corrected data to ensure order matches
corrected = integrated_obj[["integrated"]]$data %>% as.matrix() %>% t() %>% as.data.frame()
colnames(corrected) <- gsub("-", "_", colnames(corrected))
meta_corrected = integrated_obj@meta.data[,-c(1:3)]

# Write out results
corrected_df = cbind(meta_corrected, corrected)
write_parquet(corrected_df, output_file)
