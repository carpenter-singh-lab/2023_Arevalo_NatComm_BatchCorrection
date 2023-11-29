configfile: "./inputs/conf/scenario_1.json"


# Load rules
include: "rules/common.smk"
include: "rules/processing.smk"
include: "rules/map.smk"
include: "rules/sphering.smk"
include: "rules/correct.smk"
include: "rules/scib.smk"
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

REF_TYPE = ['negcon', 'nonrep']
map_pattern = f"outputs/{scenario}/metrics/{criteria}/{{workflow}}_{{method}}_map_{{reftype}}.parquet"
map_baseline_pattern = f"outputs/{scenario}/metrics/{criteria}/{{workflow}}_map_{{reftype}}.parquet"
scib_baseline_pattern = f"outputs/{scenario}/metrics/{criteria}/{{workflow}}_scib.parquet"
scib_pattern = f"outputs/{scenario}/metrics/{criteria}/{{workflow}}_{{method}}_scib.parquet"
umap_pattern = f"outputs/{scenario}/projection/{{workflow}}_{{method}}_umap.parquet"

rule all:
    input:
        expand(map_pattern, workflow=WORKFLOWS, reftype=REF_TYPE, method=METHODS),
        expand(scib_pattern, workflow=WORKFLOWS, method=METHODS),
        expand(map_baseline_pattern, workflow=WORKFLOWS, reftype=REF_TYPE),
        expand(scib_baseline_pattern, workflow=WORKFLOWS),
