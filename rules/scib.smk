METRICS = ["nmi", "ari", "asw", "silhouette_batch", "pcr_batch", "il_asw"]


rule scib_all:
    input:
        expand(
            "outputs/scenario_1/scib/mad_featselect_sphering_harmony_{metric}.bin",
            metric=METRICS,
        ),


rule scib_clustering:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/scib/{pipeline}_clusters.h5ad",
    params:
        label_key=config["label_key"],
    run:
        scib_metrics.cluster(*input, *params, *output)


rule scib_nmi:
    input:
        "outputs/{prefix}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/scib/{pipeline}_nmi.bin",
    params:
        label_key=config["label_key"],
    run:
        scib_metrics.nmi(*input, *params, *output)


rule scib_ari:
    input:
        "outputs/{prefix}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/scib/{pipeline}_ari.bin",
    params:
        label_key=config["label_key"],
    run:
        scib_metrics.ari(*input, *params, *output)


rule scib_asw:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/scib/{pipeline}_asw.bin",
    params:
        label_key=config["label_key"],
    run:
        scib_metrics.asw(*input, *params, *output)


rule scib_silhouette_batch:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/scib/{pipeline}_silhouette_batch.bin",
    params:
        label_key=config["label_key"],
        batch_key=config["batch_key"],
    run:
        scib_metrics.silhouette_batch(*input, *params, *output)


rule scib_pcr_batch:
    input:
        pre_parquet_path="outputs/{prefix}/mad.parquet",
        post_parquet_path="outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/scib/{pipeline}_pcr_batch.bin",
    params:
        batch_key=config["batch_key"],
    run:
        scib_metrics.pcr_batch(*input, *params, *output)


rule scib_il_asw:
    input:
        "outputs/{prefix}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/scib/{pipeline}_il_asw.bin",
    params:
        label_key=config["label_key"],
        batch_key=config["batch_key"],
    run:
        scib_metrics.isolated_labels_asw(*input, *params, *output)
