#!/usr/bin/env Rscript

# Suppress package startup messages for cleaner logs
suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(batchelor)
  library(SingleCellExperiment)
})

option_list = list(
  make_option(c("--input_data"), type = "character", help = "Path to input data file (parquet format)"),
  make_option(c("--batch_key"), type = "character", help = "Column name indicating batch information"),
  make_option(c("--output_path"), type = "character", help = "Path to save the corrected output file (parquet format)")
)

parser = OptionParser(option_list = option_list)
opt = parse_args(parser)

if (is.null(opt$input_data) || is.null(opt$batch_key) || is.null(opt$output_path)) {
  print_help(parser)
  stop("All three arguments (--input_data, --batch_key, --output_path) must be supplied.", call. = FALSE)
}

input_file = opt$input_data
batch_col = opt$batch_key
output_file = opt$output_path

parquet_data = read_parquet(input_file)
col_names = names(parquet_data)
metadata_cols = col_names[grepl("^Metadata_", col_names)]
features_cols = col_names[!grepl("^Metadata_", col_names)]
features = parquet_data[, features_cols]
metadata = parquet_data[, metadata_cols]
batch_info = parquet_data[, batch_col]
corrected = fastMNN(t(features), batch=batch_info)
corrected_df = cbind(metadata, reducedDim(corrected))
write_parquet(corrected_df, output_file)