'''
Functions to apply normalization methods
'''
import pandas as pd
from sklearn.base import TransformerMixin
from utils import find_meta_cols, find_feat_cols
from zca import ZCA, ZCA_corr


def apply_norm(normalizer: TransformerMixin, dframe: pd.DataFrame, field: str,
               values: list[str]):
    '''Apply normalization considering dframe[field=='value'] as training set.
    Returns the dframe normalized
    '''

    feat_cols = find_feat_cols(dframe.columns)
    meta_cols = find_meta_cols(dframe.columns)
    meta = dframe[meta_cols]
    train = dframe.loc[dframe[field].isin(values), feat_cols]
    if len(train) == 0:
        raise ValueError(f'{values} not found in {field}')
    normalizer.fit(train)

    norm_feats = normalizer.transform(dframe[feat_cols])
    norm_feats = pd.DataFrame(norm_feats,
                              index=dframe.index,
                              columns=feat_cols)
    dframe = pd.concat([meta, norm_feats], axis=1)
    return dframe


def sphere(dframe: pd.DataFrame, regularization: float, mode: str,
           column_norm: str, values_norm: list[str]):
    '''
    Apply sphering using 'negcon' as training set
    '''
    if mode == 'corr':
        normalizer = ZCA_corr(regularization=regularization)
    else:
        normalizer = ZCA(regularization=regularization)
    dframe = apply_norm(normalizer, dframe, column_norm, values_norm)
    return dframe, normalizer
