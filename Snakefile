configfile: "./inputs/conf/scenario_1.json"

wildcard_constraints:
    criteria=r"target2|prod",
    scenario=r"scenario_\d",
    pipeline=r"[_a-zA-Z.~0-9\-]*"

# Init config
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import correct
import scib_metrics
from correct import sphering
import quality_control as qc

scenario = config["scenario"]
if "COMPOUND" in config["plate_types"]:
    criteria = "prod"
else:
    criteria = "target2"


# Load rules
include: "rules/common.smk"
include: "rules/processing.smk"
include: "rules/map.smk"
include: "rules/sphering.smk"
include: "rules/correct.smk"


rule all:
    input:
        f"outputs/{scenario}/metrics/{criteria}/mad_map.parquet",
        f"outputs/{scenario}/metrics/{criteria}/mad_clip_map.parquet",
        f"outputs/{scenario}/metrics/{criteria}/mad_drop_map.parquet",
        f"outputs/{scenario}/metrics/{criteria}/mad_featselect_map.parquet",
        # f"outputs/{scenario}/metrics/{criteria}/map_mad_imputeknn.parquet",
        # f"outputs/{scenario}/metrics/{criteria}/map_mad_imputemedian.parquet",
        f"outputs/{scenario}/metrics/{criteria}/mad_int_map.parquet",
        f"outputs/{scenario}/metrics/{criteria}/mad_drop_int_map.parquet",
        f"outputs/{scenario}/metrics/{criteria}/mad_int_featselect_map.parquet",
        f"outputs/{scenario}/metrics/{criteria}/mad_drop_int_featselect_map.parquet",
        #f"outputs/{scenario}/metrics/{criteria}/mad_int_featselect_sphering_map.parquet",
        #f"outputs/{scenario}/metrics/{criteria}/mad_featselect_sphering_map.parquet",
        #f"outputs/{scenario}/metrics/{criteria}/mad_featselect_sphering_harmony_map.parquet",
        f"outputs/{scenario}/metrics/{criteria}/mad_featselect_sphering_map.parquet",
        #f"outputs/{scenario}/scib/mad_featselect_sphering_harmony_clusters.h5ad"
