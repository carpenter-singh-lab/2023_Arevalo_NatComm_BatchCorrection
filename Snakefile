configfile: "snake_params.json"


import quality_control as qc

# Load rules
include: 'rules/common.smk'
include: 'rules/mad.smk'
include: 'rules/mad_clip.smk'
include: 'rules/mad_drop.smk'
include: 'rules/mad_int.smk'
include: 'rules/mad_drop_int.smk'

rule all:
    input:
        expand(
            "outputs/{scenario}/map_target2_mad.parquet",
            scenario=config['target2_scenarios'],
        ),
        expand(
            "outputs/{scenario}/map_prod_mad.parquet",
            scenario=config['prod_scenarios'],
        ),
        expand(
            "outputs/{scenario}/map_target2_mad_drop_int.parquet",
            scenario=config['target2_scenarios'],
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_drop_int.parquet",
            scenario=config['prod_scenarios'],
        ),
