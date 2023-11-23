'''Correction methods'''
import shutil

import pandas as pd
import numpy as np

from quality_control.io import merge_parquet, split_parquet
from zca import ZCA, ZCA_corr


def log_uniform_sampling(min_=-5, max_=1, size=25, seed=[6, 12, 2022]):
    rng = np.random.default_rng(seed)
    return 10.**rng.uniform(min_, max_, size=size)


def sphering(dframe_path, mode, lambda_, column_norm, values_norm,
             sphered_path, spherer_path):
    if mode == 'corr':
        spherer = ZCA_corr(regularization=lambda_)
    elif mode == 'cov':
        spherer = ZCA(regularization=lambda_)
    else:
        raise ValueError(f'mode should be "corr" or "cov"')

    meta, vals, features = split_parquet(dframe_path)
    train_ix = meta[column_norm].isin(values_norm).values
    spherer.fit(vals[train_ix])
    vals = spherer.transform(vals).astype(np.float32)
    merge_parquet(meta, vals, features, sphered_path)
    np.savez_compressed(spherer_path, spherer=spherer)


def select_best(map_files, parquet_files, parquet_path, map_path):
    scores = []
    for map_file, parquet_file in zip(map_files, parquet_files):
        score = pd.read_parquet(map_file)['mean_average_precision'].mean()
        scores.append({
            'parquet_file': parquet_file,
            'map_file': map_file,
            'score': score
        })

    scores = pd.DataFrame(scores)
    best = scores.sort_values(by='score').iloc[-1]
    shutil.copy(best.parquet_file, parquet_path)
    shutil.copy(best.map_file, map_path)
