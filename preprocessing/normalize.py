import numpy as np
import pandas as pd
from preprocessing.io import merge_parquet, split_parquet


def mad(variant_feats_path, neg_stats_path, normalized_path):
    meta, vals, features = split_parquet(variant_feats_path)
    neg_stats = pd.read_parquet(neg_stats_path)
    neg_stats = neg_stats.query('feature in @features')

    # get counts and sort by plate
    plates, counts = np.unique(meta['Metadata_Plate'], return_counts=True)
    ix = np.argsort(meta['Metadata_Plate'])
    meta = meta.iloc[ix]
    vals = vals[ix]

    # get mad and median matrices for MAD normalization
    mads = neg_stats.pivot(index='Metadata_Plate',
                           columns='feature',
                           values='mad')
    mads = mads.loc[plates, features].values
    medians = neg_stats.pivot(index='Metadata_Plate',
                              columns='feature',
                              values='median')
    medians = medians.loc[plates, features].values

    # Get normalized features (epsilon = 0) for all plates that have MAD stats
    # -= and /= are inplace operations. i.e save memory
    vals -= np.repeat(medians, counts, axis=0)
    vals /= np.repeat(mads, counts, axis=0)

    merge_parquet(meta, vals, features, normalized_path)
