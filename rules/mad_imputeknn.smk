rule mad_imputeknn:
    input:
        "outputs/{scenario}/mad.parquet",
        "outputs/{scenario}/outliers.parquet",
    output:
        "outputs/{scenario}/mad_imputeknn.parquet",
    params:
        clip_value=config["clip_value"],
    run:
        qc.outliers.impute_knn(*input, *output)
