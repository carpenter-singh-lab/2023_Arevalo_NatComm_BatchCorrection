configfile: "./inputs/conf/scenario_2.json"


target2_scenarios = [
    "scenario_2",
    # "scenario_3",
    # "scenario_4",
    # "scenario_5",
    # "scenario_7",
]
# prod_scenarios = ["scenario_3", "scenario_5", "scenario_7"]
prod_scenarios = []


import sphering
import quality_control as qc


# Load rules
include: "rules/common.smk"
include: "rules/mad.smk"
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
        expand(
            "outputs/{scenario}/map_target2_mad.parquet",
            scenario=target2_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_prod_mad.parquet",
            scenario=prod_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_target2_mad_drop.parquet",
            scenario=target2_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_drop.parquet",
            scenario=prod_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_target2_mad_clip.parquet",
            scenario=target2_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_clip.parquet",
            scenario=prod_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_target2_mad_featselect.parquet",
            scenario=target2_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_featselect.parquet",
            scenario=prod_scenarios,
        ),
        #expand(
        #"outputs/{scenario}/map_target2_mad_imputeknn.parquet",
        #scenario=target2_scenarios,
        #),
        #expand(
        #"outputs/{scenario}/map_prod_mad_imputeknn.parquet",
        #scenario=prod_scenarios,
        #),
        expand(
            "outputs/{scenario}/map_target2_mad_int.parquet",
            scenario=target2_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_int.parquet",
            scenario=prod_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_target2_mad_int_featselect.parquet",
            scenario=target2_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_int_featselect.parquet",
            scenario=prod_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_target2_mad_drop_int.parquet",
            scenario=target2_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_drop_int.parquet",
            scenario=prod_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_target2_mad_drop_int_featselect.parquet",
            scenario=target2_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_drop_int_featselect.parquet",
            scenario=prod_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_target2_mad_int_featselect_sphering.parquet",
            scenario=target2_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_target2_mad_featselect_sphering.parquet",
            scenario=target2_scenarios,
        ),
