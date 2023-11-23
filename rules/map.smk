rule average_precision:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/ap_{criteria}_{pipeline}.parquet",
    params:
        plate_types = lambda wc: ["TARGET2"] if wc.criteria == "target2" else ["COMPOUND"]
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule mean_average_precision:
    input:
        "outputs/{prefix}/ap_{criteria}_{pipeline}.parquet",
    output:
        "outputs/{prefix}/map_{criteria}_{pipeline}.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)
