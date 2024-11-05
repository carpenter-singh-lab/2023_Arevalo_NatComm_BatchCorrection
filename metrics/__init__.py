import pandas as pd

from . import scib
from . import consistency
from .map import (
    average_precision_negcon,
    average_precision_nonrep,
    mean_average_precision,
)

DIMENSION_MAP = {
    'silhouette_batch': 'batch',
    'pcr_batch': 'batch',
    'pcr': 'batch',
    'graph_conn': 'batch',
    'kbet': 'batch',
    'lisi_batch': 'batch',
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


def concat(scib_path, negcon_path, nonrep_path, output_path):
    print(scib_path)
    scib_scores = pd.read_parquet(scib_path)

    def map_summary(path, name):
        map_scores = pd.read_parquet(path)
        map_scores.dropna(inplace=True)
        frac_p = map_scores['below_p'].sum() / len(map_scores)
        frac_q = map_scores['below_corrected_p'].sum() / len(map_scores)
        mean_map = map_scores['mean_average_precision'].mean()
        scib_scores[f'{name}_mean_map'] = mean_map
        scib_scores[f'{name}_fraction_below_p'] = frac_p
        scib_scores[f'{name}_fraction_below_corrected_p'] = frac_q

    map_summary(negcon_path, 'negcon')
    map_summary(nonrep_path, 'nonrep')

    scib_scores = scib_scores.reset_index()
    scib_scores['dimension'] = scib_scores['metric'].map(DIMENSION_MAP)
    print(scib_scores)
    scib_scores.to_parquet(output_path)
