import stats
scenarios = glob_wildcards("inputs/conf/{scenario}.json").scenario

rule all:
    input: expand("outputs/{scenario}/norm_feat_stats.parquet", scenario=scenarios)

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
