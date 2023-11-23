configfile: "./inputs/conf/scenario_2.json"

wildcard_constraints:
    criteria=r"target2|prod",
    scenario=r"scenario_\d"

import sphering
import quality_control as qc

scenario = config["scenario"]
if "COMPOUND" in config["plate_types"]:
    criteria = "prod"
else:
    criteria = "target2"


# Load rules
include: "rules/common.smk"
include: "rules/map.smk"
include: "rules/mad_clip.smk"
include: "rules/mad_drop.smk"
include: "rules/mad_imputemedian.smk"
include: "rules/mad_imputeknn.smk"
include: "rules/mad_int.smk"
include: "rules/mad_featselect.smk"
include: "rules/mad_int_featselect.smk"
include: "rules/mad_drop_int.smk"
include: "rules/mad_drop_int_featselect.smk"
include: "rules/mad_int_featselect_sphering.smk"
include: "rules/mad_featselect_sphering.smk"


rule all:
    input:
        f"outputs/{scenario}/map_{criteria}_mad.parquet",
        f"outputs/{scenario}/map_{criteria}_mad_clip.parquet",
        f"outputs/{scenario}/map_{criteria}_mad_drop.parquet",
        f"outputs/{scenario}/map_{criteria}_mad_featselect.parquet",
        f"outputs/{scenario}/map_{criteria}_mad_imputeknn.parquet",
        f"outputs/{scenario}/map_{criteria}_mad_imputemedian.parquet",
        f"outputs/{scenario}/map_{criteria}_mad_int.parquet",
        f"outputs/{scenario}/map_{criteria}_mad_drop_int.parquet",
        f"outputs/{scenario}/map_{criteria}_mad_int_featselect.parquet",
        f"outputs/{scenario}/map_{criteria}_mad_drop_int_featselect.parquet",
        f"outputs/{scenario}/map_{criteria}_mad_int_featselect_sphering.parquet",
        f"outputs/{scenario}/map_{criteria}_mad_featselect_sphering.parquet"
