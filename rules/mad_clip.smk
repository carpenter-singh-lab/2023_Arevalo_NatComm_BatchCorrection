rule mad_clip:
    input:
        "outputs/{scenario}/mad.parquet",
        "outputs/{scenario}/outliers.parquet",
    output:
        "outputs/{scenario}/mad_clip.parquet",
    params:
        clip_value=config["clip_value"],
    run:
        qc.outliers.clip_cols(*input, *params, *output)
