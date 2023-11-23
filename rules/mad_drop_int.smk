rule mad_drop_int:
    input:
        "outputs/{scenario}/mad_drop.parquet",
    output:
        "outputs/{scenario}/mad_drop_int.parquet",
    run:
        qc.transform.rank_int(*input, *output)
