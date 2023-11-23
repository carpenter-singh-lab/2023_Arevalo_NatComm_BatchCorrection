rule mad_drop:
    input:
        "outputs/{scenario}/mad.parquet",
        "outputs/{scenario}/outliers.parquet",
    output:
        "outputs/{scenario}/mad_drop.parquet",
    run:
        qc.outliers.drop_cols(*input, *output)
