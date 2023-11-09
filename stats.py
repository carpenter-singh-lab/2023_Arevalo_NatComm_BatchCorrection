from functools import partial
from itertools import chain
import json
import logging

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import median_abs_deviation
from tqdm.contrib.concurrent import thread_map

from loader import build_path, load_metadata
from utils import find_feat_cols, find_meta_cols

logging.basicConfig(format='%(levelname)s:%(asctime)s:%(name)s:%(message)s',
                    level=logging.INFO)
logger = logging.getLogger(__name__)

DMSO = 'JCP2022_033924'


def get_rows(path):
    '''Count the number of rows in a parquet file'''
    with pq.ParquetFile(path) as file:
        return file.metadata.num_rows


def prealloc_params(config):
    '''
    Get a list of paths to the parquet files and the corresponding slices
    for further concatenation
    '''
    meta = load_metadata(config['sources'], config['plate_types'])
    paths = (meta[['Metadata_Source', 'Metadata_Batch',
                   'Metadata_Plate']].drop_duplicates().apply(build_path,
                                                              axis=1)).values

    counts = thread_map(get_rows, paths)
    total = sum(counts)
    slices = np.zeros((len(paths), 2), dtype=int)
    slices[:, 1] = np.cumsum(counts)
    slices[1:, 0] = slices[:-1, 1]
    return paths, slices


def load_data(config):
    '''Load all plates given the config file'''
    paths, slices = prealloc_params(config)
    total = slices[-1, 1]

    with pq.ParquetFile(paths[0]) as f:
        meta_cols = find_meta_cols(f.schema.names)
        feat_cols = find_feat_cols(f.schema.names)
    meta = np.empty([total, len(meta_cols)], dtype='|S128')
    feats = np.empty([total, len(feat_cols)], dtype=np.float32)

    def read_parquet(params):
        path, start, end = params
        df = pd.read_parquet(path)
        meta[start:end] = df[meta_cols].values
        feats[start:end] = df[feat_cols].values

    params = np.concatenate([paths[:, None], slices], axis=1)
    thread_map(read_parquet, params)

    meta = pd.DataFrame(columns=meta_cols,
                        data=meta.astype(str),
                        dtype='category')
    dframe = pd.DataFrame(columns=feat_cols, data=feats)
    for i, col in enumerate(meta_cols[::-1]):
        dframe.insert(loc=0, column=col, value=meta[col].values)
    return dframe


def write_parquet(config_path, output_file):
    '''Write the parquet dataset given the config'''
    with open(config_path) as fread:
        config = json.load(fread)
    dframe = load_data(config)
    # Efficient merge
    meta = load_metadata(config['sources'], config['plate_types'])
    meta = (dframe[[
        'Metadata_Source', 'Metadata_Plate', 'Metadata_Well'
    ]].merge(meta,
             on=['Metadata_Source', 'Metadata_Plate', 'Metadata_Well'],
             how='left'))
    for c in meta:
        dframe[c] = meta[c].values
    dframe.to_parquet(output_file)


def get_stats(dframe: pd.DataFrame):
    feat_cols = find_feat_cols(dframe)
    dframe = dframe[feat_cols + ['Metadata_Plate']]
    median = dframe.groupby('Metadata_Plate').median().round(6)
    max_ = dframe.groupby('Metadata_Plate').max()
    min_ = dframe.groupby('Metadata_Plate').min()
    count = dframe.groupby('Metadata_Plate').count()

    mad_fn = partial(median_abs_deviation, nan_policy="omit", axis=0)
    mad = dframe.groupby('Metadata_Plate').apply(mad_fn)
    mad = pd.DataFrame(index=mad.index,
                       data=np.stack(mad.values),
                       columns=feat_cols).round(6)

    median['stat'] = 'median'
    mad['stat'] = 'mad'
    min_['stat'] = 'min'
    max_['stat'] = 'max'
    count['stat'] = 'count'

    stats = pd.concat([median, mad, min_, max_, count])
    stats.reset_index(inplace=True)
    stats = stats.melt(id_vars=['Metadata_Plate', 'stat'], var_name='feature')
    stats = stats.pivot(index=['Metadata_Plate', 'feature'],
                        columns='stat',
                        values='value')
    stats.reset_index(inplace=True)
    stats['abs_coef_var'] = ((stats['mad'] /
                              stats['median']).fillna(0).abs().replace(
                                  np.inf, 0))
    return stats


def add_metadata(stats: pd.DataFrame, meta: pd.DataFrame):
    source_map = meta[['Metadata_Source', 'Metadata_Plate']].drop_duplicates()
    source_map = source_map.set_index('Metadata_Plate').Metadata_Source
    stats['Metadata_Source'] = stats['Metadata_Plate'].map(source_map)
    stats['compartment'] = stats['feature'].apply(lambda x: x.split('_')[0])
    stats['family'] = stats['feature'].apply(
        lambda x: '_'.join(x.split('_')[:3]))


def remove_nan_infs_columns(dframe: pd.DataFrame) -> pd.DataFrame:
    feat_cols = find_feat_cols(dframe)
    withnan = dframe[feat_cols].isna().sum()[lambda x: x > 0]
    withinf = (dframe[feat_cols] == np.inf).sum()[lambda x: x > 0]
    withninf = (dframe[feat_cols] == -np.inf).sum()[lambda x: x > 0]
    redlist = set(chain(withinf.index, withnan.index, withninf.index))
    return dframe[[c for c in dframe.columns if c not in redlist]]


def write_negcon_stats(parquet_path, neg_stats_path):
    '''create statistics of negative controls for all the plates'''
    logger.info('Loading data')
    dframe = pd.read_parquet(parquet_path)
    logger.info('Removing nan and inf columns')
    dframe = remove_nan_infs_columns(dframe)
    negcon = dframe.query('Metadata_JCP2022 == @DMSO')
    logger.info('computing stats for negcons')
    neg_stats = get_stats(negcon)
    logger.info('stats done.')
    add_metadata(neg_stats, dframe[find_meta_cols(dframe)])
    neg_stats.to_parquet(neg_stats_path)


def write_variant_features(neg_stats_path, variant_feats_path):
    '''select features that have mad != 0 and abs_coef_var>1e-3 in every plate for negative controls.'''
    neg_stats = pd.read_parquet(neg_stats_path)
    neg_stats = neg_stats.query('mad!=0 and abs_coef_var>1e-3')
    variant_features = set.intersection(
        *neg_stats.groupby('Metadata_Plate')['feature'].agg(set).tolist())
    variant_features = pd.DataFrame(
        {'variant_features': list(variant_features)})
    variant_features.to_parquet(variant_feats_path)


def write_normalize_features(neg_stats_path, variant_feats_path, parquet_path,
                             norm_parquet_path):
    neg_stats = pd.read_parquet(neg_stats_path)
    variant_features = pd.read_parquet(variant_feats_path).variant_features
    neg_stats = neg_stats.query('feature in @variant_features')
    dframe = pd.read_parquet(parquet_path)

    # Compute params for MAD normalization
    mads = neg_stats.pivot(columns='feature',
                           index='Metadata_Plate',
                           values='mad')
    medians = neg_stats.pivot(columns='feature',
                              index='Metadata_Plate',
                              values='median')

    # Get normalized features with epsilon = 0 for all plates that have MAD stats
    feats = dframe.query('Metadata_Plate in @mads.index')
    fnorm = (feats.set_index('Metadata_Plate')[mads.columns] - medians) / mads
    fnorm.reset_index(drop=True, inplace=True)
    for c in find_meta_cols(feats):
        fnorm[c] = feats[c].values
    fnorm.to_parquet(norm_parquet_path)
