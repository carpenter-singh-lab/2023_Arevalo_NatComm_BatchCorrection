import numpy as np
import scipy.stats as ss
from tqdm.contrib.concurrent import thread_map

from quality_control.io import merge_parquet, split_parquet


def rank_to_normal(rank, c, n):
    '''
    Standard quantile function
    '''
    x = (rank - c) / (n - 2 * c + 1)
    return ss.norm.ppf(x)


def rank_INT(array: np.ndarray,
             c: float = 3.0 / 8,
             stochastic: bool = True,
             seed: int = 0):
    '''
    Perform rank-based inverse normal transformation in a 1d numpy array. If
    stochastic is True ties are given rank randomly, otherwise ties will share
    the same value.

    Adapted from: https://github.com/edm1/rank-based-INT/blob/85cb37bb8e0d9e71bb9e8f801fd7369995b8aee7/rank_based_inverse_normal_transformation.py
    '''
    rng = np.random.default_rng(seed=seed)

    if stochastic:
        # Shuffle
        ix = rng.permutation(len(array))
        rev_ix = np.argsort(ix)
        array = array[ix]
        # Get rank, ties are determined by their position(hence why we shuffle)
        rank = ss.rankdata(array, method="ordinal")
        rank = rank[rev_ix]
    else:
        # Get rank, ties are averaged
        rank = ss.rankdata(array, method="average")

    x = (rank - c) / (len(rank) - 2 * c + 1)
    return ss.norm.ppf(x)


def write_rank_int(normalized_path, rank_int_path):
    meta, vals, features = split_parquet(normalized_path)

    def to_normal(i):
        vals[:, i] = rank_INT(vals[:, i]).astype(np.float32)

    thread_map(to_normal, range(len(features)))
    merge_parquet(meta, vals, features, rank_int_path)
