rule mad_int_featselect:
    input:
        "outputs/{scenario}/mad_int.parquet",
    output:
        "outputs/{scenario}/mad_int_featselect.parquet",
    run:
        qc.select_features(*input, *output)
