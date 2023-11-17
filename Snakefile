configfile: "snake_params.json"


import quality_control as qc

target2_scenarios = [
    "scenario_2",
    "scenario_3",
    "scenario_4",
    "scenario_5",
    "scenario_7",
]
prod_scenarios = ["scenario_3", "scenario_5", "scenario_7"]


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
            "outputs/{scenario}/map_target2_mad_drop_int.parquet",
            scenario=target2_scenarios,
        ),
        expand(
            "outputs/{scenario}/map_prod_mad_drop_int.parquet",
            scenario=prod_scenarios,
        ),


rule write_parquet:
    input:
        "inputs/conf/{scenario}.json",
    output:
        "outputs/{scenario}/raw.parquet",
    run:
        qc.io.write_parquet(*input, *output)


rule compute_negcon_stats:
    input:
        "outputs/{scenario}/raw.parquet",
    output:
        "outputs/{scenario}/neg_stats.parquet",
    run:
        qc.stats.compute_negcon_stats(*input, *output)


rule select_variant_feats:
    input:
        "outputs/{scenario}/raw.parquet",
        "outputs/{scenario}/neg_stats.parquet",
    output:
        "outputs/{scenario}/variant_feats.parquet",
    run:
        qc.stats.select_variant_features(*input, *output)


rule mad_normalize:
    input:
        "outputs/{scenario}/variant_feats.parquet",
        "outputs/{scenario}/neg_stats.parquet",
    output:
        "outputs/{scenario}/mad.parquet",
    run:
        qc.normalize.mad(*input, *output)


rule ap_target2_mad:
    input:
        "outputs/{scenario}/mad.parquet",
    output:
        "outputs/{scenario}/ap_target2_mad.parquet",
    params:
        plate_types=["TARGET2"],
        min_replicates=2,
        max_replicates=float("inf"),
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_target2_mad:
    input:
        "outputs/{scenario}/ap_target2_mad.parquet",
    output:
        "outputs/{scenario}/map_target2_mad.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)


rule ap_prod_mad:
    input:
        "outputs/{scenario}/mad.parquet",
    output:
        "outputs/{scenario}/ap_prod_mad.parquet",
    params:
        plate_types=["COMPOUND"],
        min_replicates=2,
        max_replicates=100,  # POSCONs and DMSO have a lot more
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_prod_mad:
    input:
        "outputs/{scenario}/ap_prod_mad.parquet",
    output:
        "outputs/{scenario}/map_prod_mad.parquet",
    params:
        plate_types=["COMPOUND"],
        min_replicates=2,
        max_replicates=100,  # POSCONs and DMSO have a lot more
    run:
        qc.metrics.mean_average_precision(*input, *output)


rule rank_int:
    input:
        "outputs/{scenario}/mad.parquet",
    output:
        "outputs/{scenario}/mad_int.parquet",
    run:
        qc.transform.rank_int(*input, *output)


rule compute_norm_stats:
    input:
        "outputs/{scenario}/mad.parquet",
    output:
        "outputs/{scenario}/norm_stats.parquet",
    run:
        qc.stats.compute_stats(*input, *output)


rule iqr_outliers:
    input:
        "outputs/{scenario}/mad.parquet",
        "outputs/{scenario}/norm_stats.parquet",
    output:
        "outputs/{scenario}/outliers.parquet",
    run:
        qc.outliers.iqr(config["iqr_scale"], *input, *output)


rule clip_outliers:
    input:
        "outputs/{scenario}/mad.parquet",
        "outputs/{scenario}/outliers.parquet",
    output:
        "outputs/{scenario}/mad_clip.parquet",
    params:
        clip_value=config["clip_value"],
    run:
        qc.outliers.clip_cols(*input, *params, *output)


rule drop_outliers:
    input:
        "outputs/{scenario}/mad.parquet",
        "outputs/{scenario}/outliers.parquet",
    output:
        "outputs/{scenario}/mad_drop.parquet",
    run:
        qc.outliers.drop_cols(*input, *output)


rule drop_int:
    input:
        "outputs/{scenario}/mad_drop.parquet",
    output:
        "outputs/{scenario}/mad_drop_int.parquet",
    run:
        qc.transform.rank_int(*input, *output)


rule ap_prod_mad_drop:
    input:
        "outputs/{scenario}/mad_drop.parquet",
    output:
        "outputs/{scenario}/ap_prod_mad_drop.parquet",
    params:
        plate_types=["COMPOUND"],
        min_replicates=2,
        max_replicates=100,  # POSCONs and DMSO have a lot more
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_prod_mad_drop:
    input:
        "outputs/{scenario}/ap_prod_mad_drop.parquet",
    output:
        "outputs/{scenario}/map_prod_mad_drop.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)


rule ap_target2_mad_drop:
    input:
        "outputs/{scenario}/mad_drop.parquet",
    output:
        "outputs/{scenario}/ap_target2_mad_drop.parquet",
    params:
        plate_types=["TARGET2"],
        min_replicates=2,
        max_replicates=float("inf"),
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_target2_mad_drop:
    input:
        "outputs/{scenario}/ap_target2_mad_drop.parquet",
    output:
        "outputs/{scenario}/map_target2_mad_drop.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)


rule ap_target2_mad_clip:
    input:
        "outputs/{scenario}/mad_clip.parquet",
    output:
        "outputs/{scenario}/ap_target2_mad_clip.parquet",
    params:
        plate_types=["TARGET2"],
        min_replicates=2,
        max_replicates=float("inf"),
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_target2_mad_clip:
    input:
        "outputs/{scenario}/ap_target2_mad_clip.parquet",
    output:
        "outputs/{scenario}/map_target2_mad_clip.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)


rule ap_target2_mad_int:
    input:
        "outputs/{scenario}/mad_int.parquet",
    output:
        "outputs/{scenario}/ap_target2_mad_int.parquet",
    params:
        plate_types=["TARGET2"],
        min_replicates=2,
        max_replicates=float("inf"),
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_target2_mad_int:
    input:
        "outputs/{scenario}/ap_target2_mad_int.parquet",
    output:
        "outputs/{scenario}/map_target2_mad_int.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)


rule ap_prod_mad_int:
    input:
        "outputs/{scenario}/mad_int.parquet",
    output:
        "outputs/{scenario}/ap_prod_mad_int.parquet",
    params:
        plate_types=["COMPOUND"],
        min_replicates=2,
        max_replicates=100,  # POSCONs and DMSO have a lot more
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_prod_mad_int:
    input:
        "outputs/{scenario}/ap_prod_mad_int.parquet",
    output:
        "outputs/{scenario}/map_prod_mad_int.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)


rule ap_prod_mad_drop_int:
    input:
        "outputs/{scenario}/mad_drop_int.parquet",
    output:
        "outputs/{scenario}/ap_prod_mad_drop_int.parquet",
    params:
        plate_types=["COMPOUND"],
        min_replicates=2,
        max_replicates=100,  # POSCONs and DMSO have a lot more
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_prod_mad_drop_int:
    input:
        "outputs/{scenario}/ap_prod_mad_drop_int.parquet",
    output:
        "outputs/{scenario}/map_prod_mad_drop_int.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)


rule ap_target2_mad_drop_int:
    input:
        "outputs/{scenario}/mad_drop_int.parquet",
    output:
        "outputs/{scenario}/ap_target2_mad_drop_int.parquet",
    params:
        plate_types=["TARGET2"],
        min_replicates=2,
        max_replicates=float("inf"),
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_target2_mad_drop_int:
    input:
        "outputs/{scenario}/ap_target2_mad_drop_int.parquet",
    output:
        "outputs/{scenario}/map_target2_mad_drop_int.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)
