rule harmony:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_harmony.parquet",
    params:
        batch_key=config["batch_key"],
    run:
        correct.harmony(*input, *params, *output)


rule scanorama:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_scanorama.parquet",
    params:
        batch_key=config["batch_key"],
    run:
        correct.scanorama(*input, *params, *output)


rule pca_scanorama:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_pca_scanorama.parquet",
    params:
        batch_key=config["batch_key"],
    run:
        correct.pca_scanorama(*input, *params, *output)
