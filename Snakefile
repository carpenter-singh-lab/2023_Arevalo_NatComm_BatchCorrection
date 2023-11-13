configfile: "snake_params.json"

import stats
scenarios = glob_wildcards("inputs/conf/{scenario}.json").scenario

rule all:
    input:
        expand(f'outputs/{{scenario}}/map_target2_drop_iqr_{config["iqr_scale"]}.parquet', scenario=scenarios)

rule to_parquet:
    input: 'inputs/conf/{scenario}.json'
    output: 'outputs/{scenario}/raw.parquet'
    run: stats.write_parquet(*input, *output)

rule get_negcon_stats:
    input: 'outputs/{scenario}/raw.parquet'
    output: 'outputs/{scenario}/neg_stats.parquet'
    run: stats.write_negcon_stats(*input, *output)

rule get_variant_feats:
    input: 'outputs/{scenario}/neg_stats.parquet'
    output: 'outputs/{scenario}/variant_feats.parquet'
    run: stats.write_variant_features(*input, *output)

rule normalize_features:
    input:
        'outputs/{scenario}/raw.parquet',
        'outputs/{scenario}/neg_stats.parquet',
        'outputs/{scenario}/variant_feats.parquet'
    output: 'outputs/{scenario}/normalized.parquet'
    run:
        stats.write_normalize_features(*input, *output)

rule norm_feat_stats:
    input:
        'outputs/{scenario}/normalized.parquet',
        'outputs/{scenario}/variant_feats.parquet'
    output:
        'outputs/{scenario}/norm_feat_stats.parquet'
    run:
        stats.write_normalize_feat_stats(*input, *output)

rule iqr_outliers:
    input:
        'outputs/{scenario}/normalized.parquet',
        'outputs/{scenario}/variant_feats.parquet',
        'outputs/{scenario}/norm_feat_stats.parquet'
    output:
        'outputs/{scenario}/iqr_outliers_' f'{config["iqr_scale"]}.parquet'
    run:
        stats.write_iqr_outliers(config['iqr_scale'], *input, *output)

rule map_production_drop_iqr_outliers:
    input:
        'outputs/{scenario}/normalized.parquet',
        'outputs/{scenario}/variant_feats.parquet',
        f'outputs/{{scenario}}/iqr_outliers_{config["iqr_scale"]}.parquet'
    output:
        f'outputs/{{scenario}}/ap_prod_drop_iqr_{config["iqr_scale"]}.parquet',
        f'outputs/{{scenario}}/map_prod_drop_iqr_{config["iqr_scale"]}.parquet'
    params:
        plate_types = ['COMPOUND']
    run:
        stats.write_map_drop_outlier_cols(*input, *output, params.plate_types)


rule map_target2_drop_iqr_outliers:
    input:
        'outputs/{scenario}/normalized.parquet',
        'outputs/{scenario}/variant_feats.parquet',
        f'outputs/{{scenario}}/iqr_outliers_{config["iqr_scale"]}.parquet'
    output:
        f'outputs/{{scenario}}/ap_target2_drop_iqr_{config["iqr_scale"]}.parquet',
        f'outputs/{{scenario}}/map_target2_drop_iqr_{config["iqr_scale"]}.parquet'
    params:
        plate_types=['TARGET2'],
        min_replicates=2,
        max_replicates=float('inf')
    run:
        stats.write_map_drop_outlier_cols(*input, *output, **params)
