args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  cat("Usage: Rscript fastMNN.R input_file output_file batch_col\n")
  quit()
}

library(arrow)
library(batchelor)
library(SingleCellExperiment)
input_file <- args[1]
output_file <- args[2]
batch_col <- args[3]

parquet_data <- read_parquet(input_file)
col_names <- names(parquet_data)
metadata_cols <- col_names[grepl("^Metadata_", col_names)]
features_cols <- col_names[!grepl("^Metadata_", col_names)]
features = parquet_data[, features_cols]
metadata = parquet_data[, metadata_cols]
batch_info <- parquet_data[, batch_col]
corrected = fastMNN(t(features), batch=batch_info)
corrected_df = cbind(metadata, reducedDim(corrected))
write_parquet(corrected_df, output_file)
