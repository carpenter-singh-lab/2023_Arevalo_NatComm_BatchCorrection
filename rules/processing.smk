rule drop:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
        "outputs/{scenario}/outliers.parquet",
    output:
        "outputs/{scenario}/{pipeline}_drop.parquet",
    run:
        pp.outliers.drop_cols(*input, *output)


rule clip:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
        "outputs/{scenario}/outliers.parquet",
    output:
        "outputs/{scenario}/{pipeline}_clip.parquet",
    params:
        clip_value=config["clip_value"],
    run:
        pp.outliers.clip_cols(*input, *params, *output)


rule INT:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_int.parquet",
    run:
        pp.transform.rank_int(*input, *output)


rule featselect:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_featselect.parquet",
    run:
        pp.select_features(*input, *output)


rule imputeknn:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
        "outputs/{scenario}/outliers.parquet",
    output:
        "outputs/{scenario}/{pipeline}_imputeknn.parquet",
    params:
        clip_value=config["clip_value"],
    run:
        pp.outliers.impute_knn(*input, *output)


rule imputemedian:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
        "outputs/{scenario}/outliers.parquet",
    output:
        "outputs/{scenario}/{pipeline}_imputemedian.parquet",
    params:
        clip_value=config["clip_value"],
    run:
        pp.outliers.impute_median(*input, *output)
