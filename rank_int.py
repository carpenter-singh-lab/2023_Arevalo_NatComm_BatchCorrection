'''
Taken from: https://github.com/edm1/rank-based-INT/blob/85cb37bb8e0d9e71bb9e8f801fd7369995b8aee7/rank_based_inverse_normal_transformation.py
'''
import numpy as np
import pandas as pd
import scipy.stats as ss


def rank_INT(series, c=3.0 / 8, stochastic=True):
    rng = np.random.default_rng(seed=123)

    if stochastic:
        # Shuffle
        ix = rng.permutation(len(series))
        rev_ix = np.argsort(ix)
        series = series[ix]
        # Get rank, ties are determined by their position(hence why we shuffle)
        rank = ss.rankdata(series, method="ordinal")
        rank = rank[rev_ix]
    else:
        # Get rank, ties are averaged
        rank = ss.rankdata(series, method="average")

    x = (rank - c) / (len(rank) - 2 * c + 1)
    return ss.norm.ppf(x)


def rank_INT_pd(series, c=3.0 / 8, stochastic=True):
    """ Perform rank-based inverse normal transformation on pandas series.
        If stochastic is True ties are given rank randomly, otherwise ties will
        share the same value. NaN values are ignored.

        Args:
            param1 (pandas.Series):   Series of values to transform
            param2 (Optional[float]): Constand parameter (Bloms constant)
            param3 (Optional[bool]):  Whether to randomise rank of ties

        Returns:
            pandas.Series
    """

    # Check input
    assert (isinstance(series, pd.Series))
    assert (isinstance(c, float))
    assert (isinstance(stochastic, bool))

    # Set seed
    rng = np.random.default_rng(seed=123)

    # Take original series indexes
    orig_idx = series.index

    # Drop NaNs
    series = series.loc[~pd.isnull(series)]

    # Get ranks
    if stochastic == True:
        # Shuffle by index
        series = series.loc[rng.permutation(series.index)]
        # Get rank, ties are determined by their position in the series (hence
        # why we randomised the series)
        rank = ss.rankdata(series, method="ordinal")
    else:
        # Get rank, ties are averaged
        rank = ss.rankdata(series, method="average")

    # Convert numpy array back to series
    rank = pd.Series(rank, index=series.index)

    # Convert rank to normal distribution
    transformed = rank.apply(rank_to_normal, c=c, n=len(rank))

    return transformed[orig_idx]


def rank_to_normal(rank, c, n):
    # Standard quantile function
    x = (rank - c) / (n - 2 * c + 1)
    return ss.norm.ppf(x)
