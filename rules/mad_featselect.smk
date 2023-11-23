rule mad_featselect:
    input:
        "outputs/{scenario}/mad.parquet",
    output:
        "outputs/{scenario}/mad_featselect.parquet",
    run:
        qc.select_features(*input, *output)
