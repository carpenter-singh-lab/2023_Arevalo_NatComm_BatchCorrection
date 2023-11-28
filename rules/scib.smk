METRICS = [
    "nmi",
    "ari",
    "asw",
    "silhouette_batch",
    "pcr_batch",
    "pcr",
    "il_asw",
    "il_f1",
    "graph_conn",
    "kbet",
    "lisi_label",
    "lisi_batch",
]


rule scib_all:
    input:
        expand(
            "outputs/{{scenario}}/metrics/{{criteria}}/scib/{{pipeline}}_{metric}.bin",
            metric=METRICS,
        ),
    output:
        output_path="outputs/{scenario}/metrics/{criteria}/{pipeline}_scib.parquet",
    run:
        scib_metrics.concat(*input, **output)


rule clustering:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
    params:
        label_key=config["label_key"],
    run:
        scib_metrics.cluster(*input, *params, *output)


rule nmi:
    input:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_nmi.bin",
    params:
        label_key=config["label_key"],
    run:
        scib_metrics.nmi(*input, *params, *output)


rule ari:
    input:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_ari.bin",
    params:
        label_key=config["label_key"],
    run:
        scib_metrics.ari(*input, *params, *output)


rule asw:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_asw.bin",
    params:
        label_key=config["label_key"],
    run:
        scib_metrics.asw(*input, *params, *output)


rule silhouette_batch:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_silhouette_batch.bin",
    params:
        label_key=config["label_key"],
        batch_key=config["batch_key"],
    run:
        scib_metrics.silhouette_batch(*input, *params, *output)


rule pcr_batch:
    input:
        pre_parquet_path="outputs/{prefix}/mad.parquet",
        post_parquet_path="outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_pcr_batch.bin",
    params:
        batch_key=config["batch_key"],
    run:
        scib_metrics.pcr_batch(*input, *params, *output)


rule pcr:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_pcr.bin",
    params:
        batch_key=config["batch_key"],
    run:
        scib_metrics.pcr(*input, *params, *output)


rule il_asw:
    input:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_il_asw.bin",
    params:
        label_key=config["label_key"],
        batch_key=config["batch_key"],
    run:
        scib_metrics.isolated_labels_asw(*input, *params, *output)


rule il_f1:
    input:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_il_f1.bin",
    params:
        label_key=config["label_key"],
        batch_key=config["batch_key"],
    run:
        scib_metrics.isolated_labels_f1(*input, *params, *output)


rule graph_conn:
    input:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_graph_conn.bin",
    params:
        label_key=config["label_key"],
    run:
        scib_metrics.graph_connectivity(*input, *params, *output)


rule kbet:
    input:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_kbet.bin",
    params:
        label_key=config["label_key"],
        batch_key=config["batch_key"],
    run:
        scib_metrics.kbet(*input, *params, *output)


rule lisi_label:
    input:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_lisi_label.bin",
    params:
        label_key=config["label_key"],
    run:
        scib_metrics.lisi_label(*input, *params, *output)


rule lisi_batch:
    input:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_lisi_batch.bin",
    params:
        batch_key=config["batch_key"],
    run:
        scib_metrics.lisi_batch(*input, *params, *output)
