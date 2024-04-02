args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  cat("Usage: Rscript seurat.R input_file output_file batch_col\n")
  quit()
}
input_file <- args[1]
output_file <- args[2]
batch_col <- args[3]

library(arrow)
library(Seurat)
parquet_data <- as.data.frame(read_parquet(input_file))
col_names <- names(parquet_data)
metadata_cols <- col_names[grepl("^Metadata_", col_names)]
features_cols <- col_names[!grepl("^Metadata_", col_names)]
features = parquet_data[, features_cols]
metadata = parquet_data[, metadata_cols]
batch_info <- parquet_data[, batch_col]

batch_names <- unique(batch_info)

seurat_lists <- list()
i = 1
for (batch in split(features, batch_info)) {
  # Create a Seurat object from the batch
  raw = data.frame(batch)
  names(raw)<-gsub("_","-",names(raw))
  obj <- CreateSeuratObject(counts=t(raw))
  obj[["RNA"]]$data = t(raw)
  obj[["RNA"]]$scale.data = t(raw)
  obj = FindVariableFeatures(obj)
  obj = RunPCA(obj)
  seurat_lists[i] <- obj
  i = i + 1
}
anchor_set <- FindIntegrationAnchors(object.list = seurat_lists, reduction="rpca")
integrated_obj = IntegrateData(anchor_set)
corrected = integrated_obj[["integrated"]]$data
# TODO: Make sure order of metadata and features is the same
corrected_df = cbind(metadata, t(as.matrix(corrected)))
write_parquet(corrected_df, output_file)
