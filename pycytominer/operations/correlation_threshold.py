"""
Returns list of features such that no two features have a correlation greater than a
specified threshold
"""
import numpy as np
import pandas as pd

def get_pairwise_correlation(population_df, method="pearson"):
    """Given a population dataframe, calculate all pairwise correlations.

    Parameters
    ----------
    population_df : pandas.core.frame.DataFrame
        Includes metadata and observation features.
    method : str, default "pearson"
        Which correlation matrix to use to test cutoff.
    Returns
    -------
    list of str
        Features to exclude from the population_df.

    """

    # Get a symmetrical correlation matrix. Use numpy for non NaN/Inf matrices.
    has_nan = np.any(np.isnan(population_df.values))
    has_inf = np.any(np.isinf(population_df.values))
    if method == "pearson" and not (has_nan or has_inf):
        pop_names = population_df.columns
        data_cor_df = np.corrcoef(population_df.transpose())
        data_cor_df = pd.DataFrame(data_cor_df, index=pop_names, columns=pop_names)
    else:
        data_cor_df = population_df.corr(method=method)

    # Create a copy of the dataframe to generate upper triangle of zeros
    data_cor_natri_df = data_cor_df.copy()

    # Replace upper triangle in correlation matrix with NaN
    data_cor_natri_df = data_cor_natri_df.where(
        np.tril(np.ones(data_cor_natri_df.shape), k=-1).astype(bool)
    )

    # Acquire pairwise correlations in a long format
    # Note that we are using the NaN upper triangle DataFrame
    pairwise_df = data_cor_natri_df.stack().reset_index()
    pairwise_df.columns = ["pair_a", "pair_b", "correlation"]

    return data_cor_df, pairwise_df

def correlation_threshold(
    population_df, features="infer", samples="all", threshold=0.9, method="pearson"
):
    """Exclude features that have correlations above a certain threshold

    Parameters
    ----------
    population_df : pandas.core.frame.DataFrame
        DataFrame that includes metadata and observation features.
    features : list, default "infer"
         List of features present in the population dataframe [default: "infer"]
         if "infer", then assume cell painting features are those that start with
         "Cells_", "Nuclei_", or "Cytoplasm_".
    samples : list or str, default "all"
        List of samples to perform operation on. If "all", use all samples to calculate.
    threshold - float, default 0.9
        Must be between (0, 1) to exclude features
    method - str, default "pearson"
        indicating which correlation metric to use to test cutoff

    Returns
    -------
    excluded_features : list of str
         List of features to exclude from the population_df.
    """


    assert 0 <= threshold <= 1, "threshold variable must be between (0 and 1)"

    # Subset dataframe and calculate correlation matrix across subset features
    if samples != "all":
        population_df = population_df.loc[samples, :]

    population_df = population_df.loc[:, features]

    # Get correlation matrix and lower triangle of pairwise correlations in long format
    data_cor_df, pairwise_df = get_pairwise_correlation(
        population_df=population_df, method=method
    )

    # Get absolute sum of correlation across features
    # The lower the index, the less correlation to the full data frame
    # We want to drop features with highest correlation, so drop higher index
    variable_cor_sum = data_cor_df.abs().sum().sort_values().index

    # And subset to only variable combinations that pass the threshold
    pairwise_df = pairwise_df.query("correlation > @threshold")

    # Return an empty list if nothing is over correlation threshold
    if pairwise_df.shape[0] == 0:
        return []

    # Output the excluded features
    excluded = pairwise_df.apply(
        lambda x: determine_high_cor_pair(x, variable_cor_sum), axis="columns"
    )

    return list(set(excluded.tolist()))


def determine_high_cor_pair(correlation_row, sorted_correlation_pairs):
    """Select highest correlated variable given a correlation row with columns:
    ["pair_a", "pair_b", "correlation"]. For use in a pandas.apply().

    Parameters
    ----------
    correlation_row : pandas.core.series.series
        Pandas series of the specific feature in the pairwise_df
    sorted_correlation_pairs : pandas.DataFrame.index
        A sorted object by total correlative sum to all other features

    Returns
    -------
    The feature that has a lower total correlation sum with all other features
    """

    pair_a = correlation_row["pair_a"]
    pair_b = correlation_row["pair_b"]

    if sorted_correlation_pairs.get_loc(pair_a) > sorted_correlation_pairs.get_loc(
        pair_b
    ):
        return pair_a
    else:
        return pair_b
