import logging

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer, SimpleImputer

from quality_control.io import merge_parquet, split_parquet

logger = logging.getLogger(__name__)


def iqr(scale: float, normalized_path, stats_path, outlier_path):
    desc = pd.read_parquet(stats_path)
    meta, vals, features = split_parquet(normalized_path)

    cutoff = desc['iqr'] * scale
    lower, higher = desc['25%'] - cutoff, desc['75%'] + cutoff
    logger.info(f'Lowest/Highest threshold: {lower.min()}, {higher.min()}')
    outliers = np.logical_or(vals < lower.values, vals > higher.values)
    outliers = pd.DataFrame(outliers, columns=features)
    for c in meta:
        outliers[c] = meta[c]
    outliers.to_parquet(outlier_path)


def drop_cols(normalized_path, outlier_path, drop_outliers_path):
    '''
    Compute mAP dropping the columns with at least one outlier. It ignores DMSO
    '''
    meta, vals, features = split_parquet(normalized_path)
    mask = pd.read_parquet(outlier_path)[features].values
    no_outlier_cols = mask.sum(axis=0) == 0
    vals = vals[:, no_outlier_cols]
    features = np.asarray(features)[no_outlier_cols]
    merge_parquet(meta, vals, features, drop_outliers_path)


def clip_cols(normalized_path, outlier_path, clip_value, clip_outliers_path):
    '''
    Compute mAP clipping values to a given magnitude. It ignores DMSO
    '''
    meta, vals, features = split_parquet(normalized_path)
    mask = pd.read_parquet(outlier_path)[features].values
    vals[mask] = np.clip(vals[mask], -clip_value, clip_value)
    merge_parquet(meta, vals, features, clip_outliers_path)


def impute_median(normalized_path, outlier_path, impute_median_path):
    '''
    Impute outliers using median
    '''
    meta, vals, features = split_parquet(normalized_path)
    mask = pd.read_parquet(outlier_path)[features].values
    vals[mask] = np.nan

    imputer = SimpleImputer(copy=False, strategy='median')
    imputer.fit_transform(vals)

    merge_parquet(meta, vals, features, impute_median_path)


def impute_knn(normalized_path, outlier_path, impute_knn_path):
    '''
    Impute outliers using kNN.
    '''
    meta, vals, features = split_parquet(normalized_path)
    mask = pd.read_parquet(outlier_path)[features].values
    vals[mask] = np.nan

    imputer = KNNImputer(copy=False)
    imputer.fit_transform(vals)

    merge_parquet(meta, vals, features, impute_knn_path)
