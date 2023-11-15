configfile: "snake_params.json"

import quality_control as qc
target2_scenarios = ['scenario_2', 'scenario_3', 'scenario_4', 'scenario_5', 'scenario_7']
prod_scenarios = ['scenario_3', 'scenario_5', 'scenario_7']

rule all:
    input:
        expand('outputs/{scenario}/map_target2_drop_iqr.parquet', scenario=target2_scenarios),
        expand('outputs/{scenario}/map_prod_drop_iqr.parquet', scenario=prod_scenarios)

rule write_parquet:
    input: 'inputs/conf/{scenario}.json'
    output: 'outputs/{scenario}/raw.parquet'
    run: qc.io.write_parquet(*input, *output)

rule compute_negcon_stats:
    input: 'outputs/{scenario}/raw.parquet'
    output: 'outputs/{scenario}/neg_stats.parquet'
    run: qc.stats.compute_negcon_stats(*input, *output)

rule select_variant_feats:
    input: 'outputs/{scenario}/raw.parquet', 'outputs/{scenario}/neg_stats.parquet'
    output: 'outputs/{scenario}/variant_feats.parquet'
    run: qc.stats.select_variant_features(*input, *output)

rule mad_normalize:
    input:
        'outputs/{scenario}/variant_feats.parquet',
        'outputs/{scenario}/neg_stats.parquet',
    output: 'outputs/{scenario}/mad.parquet'
    run:
        qc.normalize.mad(*input, *output)

rule compute_norm_stats:
    input:
        'outputs/{scenario}/mad.parquet',
    output:
        'outputs/{scenario}/norm_stats.parquet'
    run:
        qc.stats.compute_stats(*input, *output)

rule rank_int:
    input:
        'outputs/{scenario}/mad.parquet',
    output:
        'outputs/{scenario}/mad_rank_int.parquet'
    run:
        qc.transform.rank_int(*input, *output)

rule iqr_outliers:
    input:
        'outputs/{scenario}/mad.parquet',
        'outputs/{scenario}/norm_stats.parquet'
    output:
        'outputs/{scenario}/iqr_outliers.parquet'
    run:
        qc.outliers.iqr(config['iqr_scale'], *input, *output)

rule clip_iqr_outliers:
    input:
        'outputs/{scenario}/mad.parquet',
        'outputs/{scenario}/iqr_outliers.parquet'
    output:
        'outputs/{scenario}/mad_clip_outlier_cols.parquet'
    params:
        clip_value = config['clip_value']
    run:
        qc.outliers.clip_cols(*input, *params, *output)

rule drop_outlier_cols:
    input:
        'outputs/{scenario}/mad.parquet',
        'outputs/{scenario}/iqr_outliers.parquet'
    output:
        'outputs/{scenario}/mad_drop_outlier_cols.parquet'
    run:
        qc.outliers.drop_cols(*input, *output)

rule ap_production_drop_iqr_outliers:
    input:
        'outputs/{scenario}/mad_drop_outlier_cols.parquet'
    output:
        'outputs/{scenario}/ap_prod_drop_iqr.parquet',
    params:
        plate_types = ['COMPOUND'],
        min_replicates = 2,
        max_replicates = 100 # POSCONs and DMSO have a lot more
    run:
        qc.evaluation.average_precision(*input, *output, **params)

rule map_production_drop_iqr_outliers:
    input:
        'outputs/{scenario}/ap_prod_drop_iqr.parquet'
    output:
        'outputs/{scenario}/map_prod_drop_iqr.parquet'
    run:
        qc.evaluation.mean_average_precision(*input, *output)

rule ap_target2_drop_iqr_outliers:
    input:
        'outputs/{scenario}/mad_drop_outlier_cols.parquet'
    output:
        'outputs/{scenario}/ap_target2_drop_iqr.parquet'
    params:
        plate_types = ['TARGET2'],
        min_replicates = 2,
        max_replicates = float('inf')
    run:
        qc.evaluation.average_precision(*input, *output, **params)

rule map_target2_drop_iqr_outliers:
    input:
        'outputs/{scenario}/ap_target2_drop_iqr.parquet'
    output:
        'outputs/{scenario}/map_target2_drop_iqr.parquet'
    run:
        qc.evaluation.mean_average_precision(*input, *output)

rule ap_target2_clip_iqr_outliers:
    input:
        'outputs/{scenario}/mad_clip_outlier_cols.parquet'
    output:
        'outputs/{scenario}/ap_target2_clip_iqr.parquet'
    params:
        plate_types = ['TARGET2'],
        min_replicates = 2,
        max_replicates = float('inf')
    run:
        qc.evaluation.average_precision(*input, *output, **params)

rule map_target2_clip_iqr_outliers:
    input:
        'outputs/{scenario}/ap_target2_clip_iqr.parquet'
    output:
        'outputs/{scenario}/map_target2_clip_iqr.parquet'
    run:
        qc.evaluation.mean_average_precision(*input, *output)

rule ap_target2_rank_int:
    input:
        'outputs/{scenario}/mad_rank_int.parquet'
    output:
        'outputs/{scenario}/ap_target2_rank_int.parquet'
    params:
        plate_types = ['TARGET2'],
        min_replicates = 2,
        max_replicates = float('inf')
    run:
        qc.evaluation.average_precision(*input, *output, **params)

rule map_target2_rank_int:
    input:
        'outputs/{scenario}/ap_target2_rank_int.parquet'
    output:
        'outputs/{scenario}/map_target2_rank_int.parquet'
    run:
        qc.evaluation.mean_average_precision(*input, *output)
