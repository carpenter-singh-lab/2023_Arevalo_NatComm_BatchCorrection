wildcard_constraints:
    criteria=r"target2|prod",
    scenario=r"scenario_\d",
    pipeline=r"[_a-zA-Z.~0-9\-]*",


# Init config
import os

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
import correct
import metrics
from correct import sphering
import preprocessing as pp
import plot

scenario = config["scenario"]
if "COMPOUND" in config["plate_types"]:
    criteria = "prod"
else:
    criteria = "target2"

rule write_parquet:
    output:
        "outputs/{scenario}/raw.parquet",
    run:
        pp.io.write_parquet(config["sources"], config["plate_types"], *output)


rule compute_negcon_stats:
    input:
        "outputs/{scenario}/raw.parquet",
    output:
        "outputs/{scenario}/neg_stats.parquet",
    run:
        pp.stats.compute_negcon_stats(*input, *output)


rule select_variant_feats:
    input:
        "outputs/{scenario}/raw.parquet",
        "outputs/{scenario}/neg_stats.parquet",
    output:
        "outputs/{scenario}/variant_feats.parquet",
    log:
        "logs/{scenario}/variant_feats.log"
    run:
        pp.stats.select_variant_features(*input, *output)


rule mad_normalize:
    input:
        "outputs/{scenario}/variant_feats.parquet",
        "outputs/{scenario}/neg_stats.parquet",
    output:
        "outputs/{scenario}/mad.parquet",
    run:
        pp.normalize.mad(*input, *output)


rule compute_norm_stats:
    input:
        "outputs/{scenario}/mad.parquet",
    output:
        "outputs/{scenario}/norm_stats.parquet",
    run:
        pp.stats.compute_stats(*input, *output)


rule iqr_outliers:
    input:
        "outputs/{scenario}/mad.parquet",
        "outputs/{scenario}/norm_stats.parquet",
    output:
        "outputs/{scenario}/outliers.parquet",
    run:
        pp.outliers.iqr(config["iqr_scale"], *input, *output)
