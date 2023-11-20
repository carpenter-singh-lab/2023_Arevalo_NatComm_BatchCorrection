rule mad_imputeknn:
    input:
        "outputs/{scenario}/mad.parquet",
        "outputs/{scenario}/outliers.parquet",
    output:
        "outputs/{scenario}/mad_imputeknn.parquet",
    params:
        clip_value=config["clip_value"],
    run:
        qc.outliers.impute_knn(*input, *output)


rule ap_target2_mad_imputeknn:
    input:
        "outputs/{scenario}/mad_imputeknn.parquet",
    output:
        "outputs/{scenario}/ap_target2_mad_imputeknn.parquet",
    params:
        plate_types=["TARGET2"],
        min_replicates=2,
        max_replicates=float("inf"),
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_target2_mad_imputeknn:
    input:
        "outputs/{scenario}/ap_target2_mad_imputeknn.parquet",
    output:
        "outputs/{scenario}/map_target2_mad_imputeknn.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)


rule ap_prod_mad_imputeknn:
    input:
        "outputs/{scenario}/mad_imputeknn.parquet",
    output:
        "outputs/{scenario}/ap_prod_mad_imputeknn.parquet",
    params:
        plate_types=["COMPOUND"],
        min_replicates=2,
        max_replicates=100,  # POSCONs and DMSO have a lot more
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_prod_mad_imputeknn:
    input:
        "outputs/{scenario}/ap_prod_mad_imputeknn.parquet",
    output:
        "outputs/{scenario}/map_prod_mad_imputeknn.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)
