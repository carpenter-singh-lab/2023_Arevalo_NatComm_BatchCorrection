rule mad_drop:
    input:
        "outputs/{scenario}/mad.parquet",
        "outputs/{scenario}/outliers.parquet",
    output:
        "outputs/{scenario}/mad_drop.parquet",
    run:
        qc.outliers.drop_cols(*input, *output)


rule ap_prod_mad_drop:
    input:
        "outputs/{scenario}/mad_drop.parquet",
    output:
        "outputs/{scenario}/ap_prod_mad_drop.parquet",
    params:
        plate_types=["COMPOUND"],
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_prod_mad_drop:
    input:
        "outputs/{scenario}/ap_prod_mad_drop.parquet",
    output:
        "outputs/{scenario}/map_prod_mad_drop.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)


rule ap_target2_mad_drop:
    input:
        "outputs/{scenario}/mad_drop.parquet",
    output:
        "outputs/{scenario}/ap_target2_mad_drop.parquet",
    params:
        plate_types=["TARGET2"],
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_target2_mad_drop:
    input:
        "outputs/{scenario}/ap_target2_mad_drop.parquet",
    output:
        "outputs/{scenario}/map_target2_mad_drop.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)
