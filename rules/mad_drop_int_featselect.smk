rule mad_drop_int_featselect:
    input:
        "outputs/{scenario}/mad_drop_int.parquet",
    output:
        "outputs/{scenario}/mad_drop_int_featselect.parquet",
    run:
        qc.select_features(*input, *output)
