'''Sphering correction'''
import shutil

import numpy as np
import pandas as pd

from preprocessing.io import merge_parquet, split_parquet
from pycytominer.operations import Spherize


def log_uniform_sampling(min_=-5, max_=3, size=25, seed=42):
    if seed is None:
        seed = [6, 12, 2022]
    rng = np.random.default_rng(seed)
    return np.round(10.**rng.uniform(min_, max_, size=size), 3)


def sphering(
    dframe_path, method, epsilon, column_norm, values_norm,
             sphered_path, spherer_path):
    spherer = Spherize(epsilon=epsilon, method=method)
    meta, vals, features = split_parquet(dframe_path)
    train_ix = meta[column_norm].isin(values_norm).values
    spherer.fit(vals[train_ix])
    vals = spherer.transform(vals).astype(np.float32)
    merge_parquet(meta, vals, features, sphered_path)
    np.savez_compressed(spherer_path, spherer=spherer)


def select_best(parquet_files,
                map_negcon_files,
                map_nonrep_files,
                ap_negcon_path,
                ap_nonrep_path,
                map_negcon_path,
                map_nonrep_path,
                parquet_path,
                ):
    scores = []
    for negcon_file, nonrep_file, parquet_file in zip(map_negcon_files,
                                                      map_nonrep_files,
                                                      parquet_files):
        negcon_score = pd.read_parquet(
            negcon_file)['mean_average_precision'].dropna().mean()
        nonrep_score = pd.read_parquet(
            nonrep_file)['mean_average_precision'].dropna().mean()
        scores.append({
            'parquet_file': parquet_file,
            'negcon_file': negcon_file,
            'nonrep_file': nonrep_file,
            'score': (negcon_score + nonrep_score) / 2
        })

    scores = pd.DataFrame(scores)
    best = scores.sort_values(by='score').iloc[-1]
    shutil.copy(best['parquet_file'], parquet_path)
    shutil.copy(best['negcon_file'], map_negcon_path)
    shutil.copy(best['nonrep_file'], map_nonrep_path)

    ap_negcon = best['negcon_file'].replace('_map_negcon.parquet','_ap_negcon.parquet')
    ap_nonrep = best['nonrep_file'].replace('_map_nonrep.parquet','_ap_nonrep.parquet')
    shutil.copy(ap_negcon, ap_negcon_path)
    shutil.copy(ap_nonrep, ap_nonrep_path)