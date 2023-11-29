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
    "sphering"
]

WORKFLOWS = [
    "mad",
    "mad_clip",
    "mad_featselect",
    "mad_int",
    "mad_int_featselect",
    "mad_drop",
    "mad_drop_featselect",
    "mad_drop",
    "mad_drop_int",
    "mad_drop_int_featselect",
]

REF_TYPE = ['negcon', 'nonrep']
map_pattern = f"outputs/{scenario}/metrics/{criteria}/{{workflow}}_{{method}}_map_{{reftype}}.parquet"
scib_pattern = f"outputs/{scenario}/metrics/{criteria}/{{workflow}}_{{method}}_scib.parquet"
mde_pattern = f"outputs/{scenario}/projection/{{workflow}}_{{method}}_mde.parquet"

rule all:
    input:
        expand(map_pattern, workflow=WORKFLOWS, reftype=REF_TYPE, method=METHODS),
        expand(scib_pattern, workflow=WORKFLOWS, method=METHODS),
