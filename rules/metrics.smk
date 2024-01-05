include: "scib.smk"
include: "map.smk"

metrics_baseline_pattern = (
    f"outputs/{scenario}/metrics/{criteria}/{{workflow}}_all_metrics.parquet"
)
metrics_pattern = (
    f"outputs/{scenario}/metrics/{criteria}/{{workflow}}_{{method}}_all_metrics.parquet"
)

rule all_metrics:
    input:
        "outputs/{scenario}/metrics/{criteria}/{pipeline}_scib.parquet",
        "outputs/{scenario}/metrics/{criteria}/{pipeline}_map_negcon.parquet",
        "outputs/{scenario}/metrics/{criteria}/{pipeline}_map_nonrep.parquet",
    output:
        "outputs/{scenario}/metrics/{criteria}/{pipeline}_all_metrics.parquet",
    run:
        metrics.concat(*input, *output)

