configfile: "./inputs/conf/scenario_1.json"


# Load rules
include: "rules/common.smk"
include: "rules/processing.smk"
include: "rules/metrics.smk"
include: "rules/correct.smk"
include: "rules/projection.smk"


METHODS = [
    "scanorama",
    "pca_scanorama",
    "mnn",
    "harmony",
    "pca_harmony",
    "combat",
]

WORKFLOWS = [
    "mad_int_featselect",
    "mad_drop_int_featselect",
]

umap_pattern = f"outputs/{scenario}/projection/{{workflow}}_{{method}}_umap.parquet"
metrics_baseline_pattern = (
    f"outputs/{scenario}/metrics/{criteria}/{{workflow}}_all_metrics.parquet"
)
metrics_pattern = (
    f"outputs/{scenario}/metrics/{criteria}/{{workflow}}_{{method}}_all_metrics.parquet"
)


rule all:
    input:
        expand(umap_pattern, workflow=WORKFLOWS, method=METHODS),
        expand(metrics_pattern, workflow=WORKFLOWS, method=METHODS),
        expand(metrics_baseline_pattern, workflow=WORKFLOWS),
