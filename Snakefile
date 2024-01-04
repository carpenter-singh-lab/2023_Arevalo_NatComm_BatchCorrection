configfile: "./inputs/conf/scenario_1.json"


# Load rules
include: "rules/common.smk"


METHODS = ["scanorama", "mnn", "harmony", "combat", "desc", "scvi"]

WORKFLOWS = [
    "mad_int_featselect",
]

PLOTS = [
    f"umap_{'batch' if scenario=='scenario_1' else 'source'}",
    "mean_all_metrics_hbarplot",
    "map_scores_barplot",
    "all_metrics_barplot",
    "cartesian",
    "results_table",
]

metrics_baseline_pattern = (
    f"outputs/{scenario}/metrics/{criteria}/{{workflow}}_all_metrics.parquet"
)
umap_baseline_pattern = f"outputs/{scenario}/projection/{{workflow}}_umap.parquet"

metrics_pattern = (
    f"outputs/{scenario}/metrics/{criteria}/{{workflow}}_{{method}}_all_metrics.parquet"
)
umap_pattern = f"outputs/{scenario}/projection/{{workflow}}_{{method}}_umap.parquet"

plots_pattern = f"outputs/{scenario}/plots/{{plot}}.svg"


include: "rules/processing.smk"
include: "rules/metrics.smk"
include: "rules/correct.smk"
include: "rules/projection.smk"
include: "rules/plots.smk"


rule all:
    input:
        expand(metrics_baseline_pattern, workflow=WORKFLOWS),
        expand(umap_baseline_pattern, workflow=WORKFLOWS),
        expand(metrics_pattern, workflow=WORKFLOWS, method=METHODS),
        expand(umap_pattern, workflow=WORKFLOWS, method=METHODS),
        expand(plots_pattern, plot=PLOTS),
