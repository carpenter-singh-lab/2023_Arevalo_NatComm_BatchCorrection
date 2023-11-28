rule average_precision_negcon:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/metrics/{criteria}/{pipeline}_ap_negcon.parquet",
    params:
        plate_types=lambda wc: ["TARGET2"] if wc.criteria == "target2" else ["COMPOUND"],
    run:
        qc.metrics.average_precision_negcon(*input, *output, **params)


rule average_precision_nonrep:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/metrics/{criteria}/{pipeline}_ap_nonrep.parquet",
    params:
        plate_types=lambda wc: ["TARGET2"] if wc.criteria == "target2" else ["COMPOUND"],
    run:
        qc.metrics.average_precision_nonrep(*input, *output, **params)

rule mean_average_precision:
    input:
        "outputs/{prefix}/metrics/{criteria}/{pipeline}_ap_{reftype}.parquet",
    output:
        "outputs/{prefix}/metrics/{criteria}/{pipeline}_map_{reftype}.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)
