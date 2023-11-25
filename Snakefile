configfile: "./inputs/conf/scenario_1.json"


# Load rules
include: "rules/common.smk"
include: "rules/processing.smk"
include: "rules/map.smk"
include: "rules/sphering.smk"
include: "rules/correct.smk"
include: "rules/scib.smk"


rule all:
    input:
        f"outputs/{scenario}/metrics/{criteria}/mad_featselect_sphering_{criteria}_harmony_map.parquet",
        f"outputs/{scenario}/metrics/mad_featselect_sphering_{criteria}_harmony_scib.parquet",
