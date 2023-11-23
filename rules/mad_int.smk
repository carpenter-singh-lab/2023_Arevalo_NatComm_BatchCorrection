rule mad_int:
    input:
        "outputs/{scenario}/mad.parquet",
    output:
        "outputs/{scenario}/mad_int.parquet",
    run:
        qc.transform.rank_int(*input, *output)
