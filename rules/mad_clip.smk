rule mad_clip:
    input:
        "outputs/{scenario}/mad.parquet",
        "outputs/{scenario}/outliers.parquet",
    output:
        "outputs/{scenario}/mad_clip.parquet",
    params:
        clip_value=config["clip_value"],
    run:
        qc.outliers.clip_cols(*input, *params, *output)


rule ap_target2_mad_clip:
    input:
        "outputs/{scenario}/mad_clip.parquet",
    output:
        "outputs/{scenario}/ap_target2_mad_clip.parquet",
    params:
        plate_types=["TARGET2"],
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_target2_mad_clip:
    input:
        "outputs/{scenario}/ap_target2_mad_clip.parquet",
    output:
        "outputs/{scenario}/map_target2_mad_clip.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)


rule ap_prod_mad_clip:
    input:
        "outputs/{scenario}/mad_clip.parquet",
    output:
        "outputs/{scenario}/ap_prod_mad_clip.parquet",
    params:
        plate_types=["COMPOUND"],
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_prod_mad_clip:
    input:
        "outputs/{scenario}/ap_prod_mad_clip.parquet",
    output:
        "outputs/{scenario}/map_prod_mad_clip.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)
