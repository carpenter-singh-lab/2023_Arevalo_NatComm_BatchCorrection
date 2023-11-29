import pandas as pd

from . import scib
from .map import (
    average_precision_negcon,
    average_precision_nonrep,
    mean_average_precision,
)


def concat(scib_path, negcon_path, nonrep_path, output_path):
    scores = pd.read_parquet(scib_path).set_index('metric')['score']

    def map_summary(path, name):
        map_scores = pd.read_parquet(path)
        map_scores.dropna(inplace=True)
        frac_p = map_scores.above_p_threshold.sum() / len(map_scores)
        frac_q = map_scores.above_p_threshold.sum() / len(map_scores)
        mean_map = map_scores['mean_average_precision'].mean()
        scores[f'{name}_mean_map'] = mean_map
        scores[f'{name}_fraction_positive_p'] = frac_p
        scores[f'{name}_fraction_positive_q'] = frac_q

    map_summary(negcon_path, 'negcon')
    map_summary(nonrep_path, 'nonrep')

    scores.reset_index().to_parquet(output_path)
