BBKNN_METRICS = ["graph_conn", "kbet", "lisi_label", "lisi_batch", "ari", "nmi"]

rule bbknn_clustering:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_bbknn_clusters.h5ad",
    log:
        "logs/{prefix}/{criteria}_{pipeline}_bbknn_clustering.log",
    params:
        batch_key=config["batch_key"],
    run:
        correct.bbknn.clustering(*input, *output, params.batch_key)


rule bbknn_all:
    input:
        expand(
            "outputs/{{scenario}}/metrics/{{criteria}}/scib/{{pipeline}}_{metric}.npy",
            metric=BBKNN_METRICS,
        ),
    output:
        output_path="outputs/{scenario}/metrics/{criteria}/{pipeline}_bbknn_scib.parquet",
    log:
        "logs/{scenario}/{criteria}_{pipeline}_bbknn_all.log",
    run:
        metrics.scib.concat(*input, **output)


ruleorder: bbknn_clustering > clustering
ruleorder: bbknn_all > scib_all
