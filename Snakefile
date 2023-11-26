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
        # "scanorama",
    "pca_scanorama",
    "mnn",
    "harmony",
    "pca_harmony",
]
map_pattern = f"outputs/{scenario}/metrics/{criteria}/mad_featselect_sphering_{criteria}_{{method}}_map.parquet"
scib_pattern = f"outputs/{scenario}/metrics/{criteria}/mad_featselect_sphering_{criteria}_{{method}}_scib.parquet"
mde_pattern = f"outputs/{scenario}/projection/mad_featselect_sphering_{criteria}_{{method}}_mde.parquet"


rule all:
    input:
        expand(map_pattern, method=METHODS),
        expand(scib_pattern, method=METHODS),
        expand(mde_pattern, method=METHODS),
