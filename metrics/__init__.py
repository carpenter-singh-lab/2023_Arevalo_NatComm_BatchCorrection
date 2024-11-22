import pandas as pd

from . import scib
from . import consistency
from .map import (
    # average_precision_negcon,
    # average_precision_nonrep,
    mean_average_precision,
)
import numpy as np

DIMENSION_MAP = {
    'silhouette_batch': 'batch',
    'pcr_batch': 'labelfree',
    'pcr': 'labelfree',
    'graph_conn': 'batch',
    'kbet': 'labelfree',
    'lisi_batch': 'labelfree',
    'lisi_label': 'bio',
    'negcon_mean_map': 'bio',
    'negcon_fraction_below_p': 'bio',
    'negcon_fraction_below_corrected_p': 'bio',
    'nonrep_mean_map': 'bio',
    'nonrep_fraction_below_p': 'bio',
    'nonrep_fraction_below_corrected_p': 'bio',
    'nmi': 'bio',
    'ari': 'bio',
    'asw': 'bio',
    'il_f1': 'bio',
    'il_asw': 'bio',
}


def concat_and_average_all_metrics(
    scib_path: str,
    map_path: str, 
    metrics_redlist: list[str],
    methods_redlist: list[str],
    output_path: str
):

    scib_scores = pd.read_parquet(scib_path)
    map_scores = pd.read_parquet(map_path)

    expected_cols = [
        "eval_key",
        "method",
        "metric",
        "value",
        "metric_type",
    ]

    assert all(col in scib_scores.columns for col in expected_cols)
    assert all(col in map_scores.columns for col in expected_cols)

    # assert correct col order
    scib_scores = scib_scores[expected_cols]
    map_scores = map_scores[expected_cols]

    # concatenate and sort scores
    all_scores = pd.concat([scib_scores, map_scores], ignore_index=True)
    all_scores = all_scores.sort_values(by=["eval_key", "method", "metric"])

    # remove metrics and methods we don't want to include
    all_scores = all_scores[~all_scores["metric"].isin(metrics_redlist)]
    all_scores = all_scores[~all_scores["method"].isin(methods_redlist)]

    # aggregate columns
    aggregate_scores = []

    for eval_key in all_scores.eval_key.unique():
        for method in all_scores.method.unique():
            per_method_means = {}
            for group in ["batch", "bio"]:
                local_scores = all_scores[
                    (all_scores.eval_key == eval_key) 
                    & (all_scores.method == method)
                    & (all_scores.metric_type == ("batch_correction" if group == "batch" else "bio_conservation"))
                ]
                per_method_means[group] = np.mean(local_scores.value)
                aggregate_scores.append({
                    "eval_key": eval_key,
                    "method": method,
                    "metric": f"mean_{group}",
                    "value": per_method_means[group],
                    "metric_type": "aggregate_score",
                })
            
            aggregate_scores.append({
                "eval_key": eval_key,
                "method": method,
                "metric": "mean_overall",
                "value": 0.4 * per_method_means["batch"] + 0.6 * per_method_means["bio"],
                "metric_type": "aggregate_score",
            })
            
    scores_with_means = pd.concat([all_scores, pd.DataFrame(aggregate_scores)])
    scores_with_means = scores_with_means.sort_values(by=["eval_key", "method", "metric"])
    scores_with_means = scores_with_means.reset_index(drop=True)
    
    scores_with_means.reset_index(drop=True).to_parquet(output_path, index=False)
        