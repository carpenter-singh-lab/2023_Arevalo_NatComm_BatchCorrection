rule mad_imputemedian:
    input:
        "outputs/{scenario}/mad.parquet",
        "outputs/{scenario}/outliers.parquet",
    output:
        "outputs/{scenario}/mad_imputemedian.parquet",
    params:
        clip_value=config["clip_value"],
    run:
        qc.outliers.impute_median(*input, *output)
