SCIB_METRICS = [
    "nmi",
    "ari",
    # "asw",
    "silhouette_batch",
    "pcr_batch",
    "pcr",
    "graph_conn",
    # "kbet",
    # "lisi_label",
    # "lisi_batch",
]


rule scib_all:
    input:
        metric_paths=
        [expand(
            "outputs/{{scenario}}/metrics/{{criteria}}/scib/{{pipeline}}_{metric}.npy",
            metric=SCIB_METRICS,
        )] + [
        expand(
            "outputs/{{scenario}}/metrics/{{criteria}}/scib/{{pipeline}}_asw_{key}.npy",
            key=config["eval_key"],
        )],
    params:
        scib_metrics=SCIB_METRICS,
        eval_keys=config["eval_key"],
    output:
        output_path="outputs/{scenario}/metrics/{criteria}/{pipeline}_scib.parquet",
    run:
        metrics.scib.concat_metrics_per_method_into_df(input.metric_paths, params.scib_metrics, params.eval_keys, output.output_path)


rule clustering:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
    run:
        metrics.scib.cluster(*input, *output)


rule nmi:
    input:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_nmi.npy",
    params:
        label_key=config["label_key"],
    run:
        metrics.scib.nmi(*input, *params, *output)


rule ari:
    input:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_ari.npy",
    params:
        label_key=config["label_key"],
    run:
        metrics.scib.ari(*input, *params, *output)


rule asw:
    input:
        parquet_path="outputs/{prefix}/{pipeline}.parquet",
    output:
        asw_path="outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_asw_{key}.npy",
    params:
        label_key="{key}",
    run:
        metrics.scib.asw(input.parquet_path, params.label_key, output.asw_path)


rule silhouette_batch:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_silhouette_batch.npy",
    params:
        label_key=config["label_key"],
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    run:
        metrics.scib.silhouette_batch(*input, *params, *output)


rule pcr_batch:
    input:
        pre_parquet_path="outputs/{prefix}/mad_int_featselect.parquet",
        post_parquet_path="outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_pcr_batch.npy",
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    run:
        metrics.scib.pcr_batch(*input, *params, *output)


rule pcr:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_pcr.npy",
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    run:
        metrics.scib.pcr(*input, *params, *output)


rule il_asw:
    input:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_il_asw.npy",
    params:
        label_key=config["label_key"],
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    run:
        metrics.scib.isolated_labels_asw(*input, *params, *output)


rule il_f1:
    input:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_il_f1.npy",
    params:
        label_key=config["label_key"],
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    run:
        metrics.scib.isolated_labels_f1(*input, *params, *output)


rule graph_conn:
    input:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
    output:
        "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_graph_conn.npy",
    params:
        label_key=config["label_key"],
    run:
        metrics.scib.graph_connectivity(*input, *params, *output)


# rule kbet:
#     input:
#         "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
#     output:
#         "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_kbet.npy",
#     params:
#         label_key=config["label_key"],
#         batch_key=config["batch_key"],
#     run:
#         metrics.scib.kbet(*input, *params, *output)


# rule lisi_label:
#     input:
#         "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
#     output:
#         "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_lisi_label.npy",
#     params:
#         label_key=config["label_key"],
#     run:
#         metrics.scib.lisi_label(*input, *params, *output)


# rule lisi_batch:
#     input:
#         "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_clusters.h5ad",
#     output:
#         "outputs/{prefix}/metrics/{criteria}/scib/{pipeline}_lisi_batch.npy",
#     params:
#         batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
#     run:
#         metrics.scib.lisi_batch(*input, *params, *output)
