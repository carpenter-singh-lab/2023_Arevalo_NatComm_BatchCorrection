configfile: "snake_params.json"


import quality_control as qc

# Load rules
include: 'rules/common.smk'
include: 'rules/mad.smk'
include: 'rules/mad_clip.smk'
include: 'rules/mad_drop.smk'
include: 'rules/mad_imputemedian.smk'
include: 'rules/mad_imputeknn.smk'
include: 'rules/mad_int.smk'
include: 'rules/mad_featselect.smk'
include: 'rules/mad_int_featselect.smk'
include: 'rules/mad_drop_int.smk'
include: 'rules/mad_drop_int_featselect.smk'

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
            "outputs/{scenario}/map_target2_mad_drop.parquet",
            scenario=config['target2_scenarios'],
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_drop.parquet",
            scenario=config['prod_scenarios'],
        ),

        expand(
            "outputs/{scenario}/map_target2_mad_clip.parquet",
            scenario=config['target2_scenarios'],
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_clip.parquet",
            scenario=config['prod_scenarios'],
        ),

        expand(
            "outputs/{scenario}/map_target2_mad_featselect.parquet",
            scenario=config['target2_scenarios'],
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_featselect.parquet",
            scenario=config['prod_scenarios'],
        ),

        expand(
            "outputs/{scenario}/map_target2_mad_imputemedian.parquet",
            scenario=config['target2_scenarios'],
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_imputemedian.parquet",
            scenario=config['prod_scenarios'],
        ),

        expand(
            "outputs/{scenario}/map_target2_mad_imputeknn.parquet",
            scenario=config['target2_scenarios'],
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_imputeknn.parquet",
            scenario=config['prod_scenarios'],
        ),

        expand(
            "outputs/{scenario}/map_target2_mad_int.parquet",
            scenario=config['target2_scenarios'],
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_int.parquet",
            scenario=config['prod_scenarios'],
        ),

        expand(
            "outputs/{scenario}/map_target2_mad_int_featselect.parquet",
            scenario=config['target2_scenarios'],
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_int_featselect.parquet",
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

        expand(
            "outputs/{scenario}/map_target2_mad_drop_int_featselect.parquet",
            scenario=config['target2_scenarios'],
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_drop_int_featselect.parquet",
            scenario=config['prod_scenarios'],
        ),
