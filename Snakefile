configfile: "./inputs/conf/scenario_1.json"


# Load rules
include: "rules/common.smk"
include: "rules/processing.smk"
include: "rules/metrics.smk"
include: "rules/sphering.smk"
include: "rules/correct.smk"
include: "rules/projection.smk"


METHODS = [
    "scanorama",
    "pca_scanorama",
    "mnn",
    "harmony",
    "pca_harmony",
]

WORKFLOWS = [
    "mad",
    "mad_clip",
    "mad_featselect",
    "mad_int",
    "mad_int_featselect",
    "mad_int_featselect_sphering",
    "mad_drop",
    "mad_drop_featselect",
    "mad_drop",
    "mad_drop_int",
    "mad_drop_int_featselect",
    "mad_drop_int_featselect_sphering",
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
        expand(metrics_pattern, workflow=WORKFLOWS, method=METHODS),
        expand(metrics_baseline_pattern, workflow=WORKFLOWS),
