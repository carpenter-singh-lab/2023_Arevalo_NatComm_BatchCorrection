rule average_precision:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/metrics/{criteria}/{pipeline}_ap.parquet",
    params:
        plate_types = lambda wc: ["TARGET2"] if wc.criteria == "target2" else ["COMPOUND"]
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule mean_average_precision:
    input:
        "outputs/{prefix}/metrics/{criteria}/{pipeline}_ap.parquet"
    output:
        "outputs/{prefix}/metrics/{criteria}/{pipeline}_map.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)
