from functools import partial
from itertools import chain
import json
import logging

from copairs.map import aggregate, run_pipeline
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import median_abs_deviation
from tqdm.contrib.concurrent import thread_map

from loader import build_path, load_metadata
from utils import find_feat_cols, find_meta_cols

logger = logging.getLogger(__name__)


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

    counts = thread_map(get_rows, paths, leave=False, desc='counts')
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

    meta = pd.DataFrame(data=meta.astype(str),
                        columns=meta_cols,
                        dtype='category')
    dframe = pd.DataFrame(columns=feat_cols, data=feats)
    for i, col in enumerate(meta_cols):
        dframe[col] = meta[col]
    return dframe


def write_parquet(config_path, output_file):
    '''Write the parquet dataset given the config'''
    with open(config_path) as fread:
        config = json.load(fread)
    dframe = load_data(config)
    # Efficient merge
    meta = load_metadata(config['sources'], config['plate_types'])
    foreign_key = ['Metadata_Source', 'Metadata_Plate', 'Metadata_Well']
    meta = dframe[foreign_key].merge(meta, on=foreign_key, how='left')
    for c in meta:
        dframe[c] = meta[c].astype('category')
    dframe.to_parquet(output_file)


def get_feat_stats(dframe: pd.DataFrame, features=None):
    '''Get statistics per each feature'''
    if features is None:
        features = find_feat_cols(dframe)
    desc = thread_map(lambda x: dframe[x].describe(), features, leave=False)
    desc = pd.DataFrame(desc)
    desc['iqr'] = desc['75%'] - desc['25%']
    return desc


def get_plate_stats(dframe: pd.DataFrame):
    mad_fn = partial(median_abs_deviation, nan_policy="omit", axis=0)

    feat_cols = find_feat_cols(dframe)
    dframe = dframe[feat_cols + ['Metadata_Plate']]
    median = dframe.groupby('Metadata_Plate', observed=True).median()
    max_ = dframe.groupby('Metadata_Plate', observed=True).max()
    min_ = dframe.groupby('Metadata_Plate', observed=True).min()
    count = dframe.groupby('Metadata_Plate', observed=True).count()
    mad = dframe.groupby('Metadata_Plate', observed=True).apply(mad_fn)
    mad = pd.DataFrame(index=mad.index,
                       data=np.stack(mad.values),
                       columns=feat_cols)

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
    stats = stats.astype({
        'min': np.float32,
        'max': np.float32,
        'count': np.int32,
        'median': np.float32,
        'mad': np.float32,
        'abs_coef_var': np.float32,
        'feature': 'category'
    })
    return stats


def add_metadata(stats: pd.DataFrame, meta: pd.DataFrame):
    source_map = meta[['Metadata_Source', 'Metadata_Plate']].drop_duplicates()
    source_map = source_map.set_index('Metadata_Plate').Metadata_Source
    stats['Metadata_Source'] = stats['Metadata_Plate'].map(source_map)
    parts = stats['feature'].str.split('_', expand=True)
    stats['compartment'] = parts[0].astype('category')
    stats['family'] = parts[range(3)].apply('_'.join,
                                            axis=1).astype('category')


def remove_nan_infs_columns(dframe: pd.DataFrame) -> pd.DataFrame:
    feat_cols = find_feat_cols(dframe)
    withnan = dframe[feat_cols].isna().sum()[lambda x: x > 0]
    withinf = (dframe[feat_cols] == np.inf).sum()[lambda x: x > 0]
    withninf = (dframe[feat_cols] == -np.inf).sum()[lambda x: x > 0]
    redlist = set(chain(withinf.index, withnan.index, withninf.index))
    return dframe[[c for c in dframe.columns if c not in redlist]]


def write_negcon_stats(parquet_path, neg_stats_path):
    '''create statistics of negative controls platewise'''
    logger.info('Loading data')
    dframe = pd.read_parquet(parquet_path)
    logger.info('Removing nan and inf columns')
    dframe = remove_nan_infs_columns(dframe)
    negcon = dframe.query('Metadata_JCP2022 == "DMSO"')
    logger.info('computing stats for negcons')
    neg_stats = get_plate_stats(negcon)
    logger.info('stats done.')
    add_metadata(neg_stats, dframe[find_meta_cols(dframe)])
    neg_stats.to_parquet(neg_stats_path)


def write_normalize_feat_stats(normalized_path, stats_path):
    dframe = pd.read_parquet(normalized_path)
    fea_stats = get_feat_stats(dframe)
    fea_stats.to_parquet(stats_path)


def write_variant_features(neg_stats_path, variant_feats_path):
    '''select features that have mad != 0 and abs_coef_var>1e-3 in every plate for negative controls.'''
    neg_stats = pd.read_parquet(neg_stats_path)
    neg_stats = neg_stats.query('mad!=0 and abs_coef_var>1e-3')
    groups = neg_stats.groupby('Metadata_Plate', observed=True)['feature']
    variant_features = set.intersection(*groups.agg(set).tolist())
    variant_features = pd.DataFrame(
        {'variant_features': sorted(variant_features)})
    variant_features.to_parquet(variant_feats_path)


def write_normalize_features(parquet_path, neg_stats_path, variant_feats_path,
                             normalized_path):
    dframe = pd.read_parquet(parquet_path)
    neg_stats = pd.read_parquet(neg_stats_path)
    variant_features = np.ravel(pd.read_parquet(variant_feats_path))

    # Use only plates with variant features
    neg_stats = neg_stats.query('feature in @variant_features')
    dframe = dframe.query('Metadata_Plate in @neg_stats.Metadata_Plate')
    dframe = dframe.sort_values(by='Metadata_Plate')
    plate_counts = dframe['Metadata_Plate'].value_counts(sort=False)
    plate_counts = plate_counts[plate_counts > 0]
    plates, counts = plate_counts.index, plate_counts.values
    meta = dframe[find_meta_cols(dframe)]
    feats = dframe[variant_features].values
    del dframe

    # get mad and median matrices for MAD normalization
    mads = neg_stats.pivot(index='Metadata_Plate',
                           columns='feature',
                           values='mad')
    mads = mads.loc[plates, variant_features].values
    medians = neg_stats.pivot(index='Metadata_Plate',
                              columns='feature',
                              values='median')
    medians = medians.loc[plates, variant_features].values

    # Get normalized features (epsilon = 0) for all plates that have MAD stats
    # -= and /= are inplace operations. i.e save memory
    feats -= np.repeat(medians, counts, axis=0)
    feats /= np.repeat(mads, counts, axis=0)

    dframe = pd.DataFrame(feats, columns=variant_features)
    for c in meta:
        dframe[c] = meta[c].values
    dframe.to_parquet(normalized_path)


def load_separated(dframe_path, features=None, return_names=False):
    dframe = pd.read_parquet(dframe_path)
    if features is None:
        features = find_feat_cols(dframe)
    vals = np.empty((len(features), len(dframe)), dtype=np.float32)
    for i, c in enumerate(features):
        vals[i] = dframe[c]
    meta = dframe[find_meta_cols(dframe)].copy()
    if return_names:
        return meta, vals, features
    return meta, vals


def write_iqr_outliers(
    scale: float,
    normalized_path,
    stats_path,
    outlier_path,
):
    desc = pd.read_parquet(stats_path)
    meta, vals, features = load_separated(normalized_path, return_names=True)

    cutoff = desc['iqr'] * scale
    lower, higher = desc['25%'] - cutoff, desc['75%'] + cutoff
    logger.info(f'Lowest/Highest threshold: {lower.min()}, {higher.min()}')
    outliers = np.logical_or(vals.T < lower.values, vals.T > higher.values)
    outliers = pd.DataFrame(outliers, columns=features)
    for c in meta:
        outliers[c] = meta[c]
    outliers.to_parquet(outlier_path)


def write_map_drop_outlier_cols(normalized_path,
                                outlier_path,
                                ap_path,
                                map_path,
                                plate_types,
                                min_replicates,
                                max_replicates,
                                ignore_dmso=True):
    '''
    Compute mAP dropping the columns with at least one outlier. It ignores DMSO
    '''
    meta, vals, features = load_separated(normalized_path, return_names=True)
    mask = pd.read_parquet(outlier_path)[features].values

    ix = meta['Metadata_PlateType'].isin(plate_types)

    # Select compounds to be used in mAP computation
    valid_cmpd = meta.loc[ix, 'Metadata_JCP2022'].value_counts()
    if ignore_dmso:
        valid_cmpd.drop('DMSO', inplace=True)
    valid_cmpd = valid_cmpd[valid_cmpd.between(min_replicates, max_replicates)]
    valid_cmpd = valid_cmpd.index

    ix = (ix & meta['Metadata_JCP2022'].isin(valid_cmpd)).values

    result = run_pipeline(
        meta[ix],
        vals[mask.sum(axis=0) == 0].T[ix],
        pos_sameby=['Metadata_JCP2022'],
        pos_diffby=[],
        neg_sameby=['Metadata_Plate'],
        neg_diffby=['Metadata_JCP2022'],
        null_size=10000,
        batch_size=20000,
        seed=0,
    )
    agg_result = aggregate(result, 'Metadata_JCP2022', threshold=0.05)
    agg_result = agg_result.query('Metadata_JCP2022 in @valid_cmpd')
    result.to_parquet(ap_path)
    agg_result.to_parquet(map_path)


def write_map_clip_outlier_cols(normalized_path,
                                outlier_path,
                                ap_path,
                                map_path,
                                plate_types,
                                min_replicates,
                                max_replicates,
                                clip_value,
                                ignore_dmso=True):
    '''
    Compute mAP clipping values to a given magnitude. It ignores DMSO
    '''
    meta, vals, features = load_separated(normalized_path, return_names=True)
    mask = pd.read_parquet(outlier_path)[features].values

    ix = meta['Metadata_PlateType'].isin(plate_types)

    # Select compounds to be used in mAP computation
    valid_cmpd = meta.loc[ix, 'Metadata_JCP2022'].value_counts()
    if ignore_dmso:
        valid_cmpd.drop('DMSO', inplace=True)
    valid_cmpd = valid_cmpd[valid_cmpd.between(min_replicates, max_replicates)]
    valid_cmpd = valid_cmpd.index

    ix = (ix & meta['Metadata_JCP2022'].isin(valid_cmpd)).values
    np.clip(vals, -clip_value, clip_value, out=vals)

    result = run_pipeline(
        meta[ix],
        vals.T[ix],
        pos_sameby=['Metadata_JCP2022'],
        pos_diffby=[],
        neg_sameby=['Metadata_Plate'],
        neg_diffby=['Metadata_JCP2022'],
        null_size=10000,
        batch_size=20000,
        seed=0,
    )
    agg_result = aggregate(result, 'Metadata_JCP2022', threshold=0.05)
    agg_result = agg_result.query('Metadata_JCP2022 in @valid_cmpd')
    result.to_parquet(ap_path)
    agg_result.to_parquet(map_path)


def write_map_impute_median_outlier_cols(normalized_path,
                                         outlier_path,
                                         ap_path,
                                         map_path,
                                         plate_types,
                                         min_replicates,
                                         max_replicates,
                                         ignore_dmso=True):
    '''
    Compute mAP clipping values to a given magnitude. It ignores DMSO
    '''
    meta, vals, features = load_separated(normalized_path, return_names=True)
    mask = pd.read_parquet(outlier_path)[features].values

    ix = meta['Metadata_PlateType'].isin(plate_types)

    # Select compounds to be used in mAP computation
    valid_cmpd = meta.loc[ix, 'Metadata_JCP2022'].value_counts()
    if ignore_dmso:
        valid_cmpd.drop('DMSO', inplace=True)
    valid_cmpd = valid_cmpd[valid_cmpd.between(min_replicates, max_replicates)]
    valid_cmpd = valid_cmpd.index

    from sklearn.impute import SimpleImputer
    imputer = SimpleImputer(copy=False, strategy='median')
    vals = vals.T[ix]
    vals[mask] = np.nan
    imputer.fit_transform(vals)


    ix = (ix & meta['Metadata_JCP2022'].isin(valid_cmpd)).values
    result = run_pipeline(
        meta[ix],
        vals[ix],
        pos_sameby=['Metadata_JCP2022'],
        pos_diffby=[],
        neg_sameby=['Metadata_Plate'],
        neg_diffby=['Metadata_JCP2022'],
        null_size=10000,
        batch_size=20000,
        seed=0,
    )
    agg_result = aggregate(result, 'Metadata_JCP2022', threshold=0.05)
    agg_result = agg_result.query('Metadata_JCP2022 in @valid_cmpd')
    result.to_parquet(ap_path)
    agg_result.to_parquet(map_path)


def write_map_impute_knn_outlier_cols(normalized_path,
                                         outlier_path,
                                         ap_path,
                                         map_path,
                                         plate_types,
                                         min_replicates,
                                         max_replicates,
                                         ignore_dmso=True):
    '''
    Compute mAP clipping values to a given magnitude. It ignores DMSO
    '''
    meta, vals, features = load_separated(normalized_path, return_names=True)
    mask = pd.read_parquet(outlier_path)[features].values

    ix = meta['Metadata_PlateType'].isin(plate_types)

    # Select compounds to be used in mAP computation
    valid_cmpd = meta.loc[ix, 'Metadata_JCP2022'].value_counts()
    if ignore_dmso:
        valid_cmpd.drop('DMSO', inplace=True)
    valid_cmpd = valid_cmpd[valid_cmpd.between(min_replicates, max_replicates)]
    valid_cmpd = valid_cmpd.index

    from sklearn.impute import KNNImputer
    imputer = KNNImputer(copy=False)
    vals = vals.T[ix]
    vals[mask] = np.nan
    imputer.fit_transform(vals)


    ix = (ix & meta['Metadata_JCP2022'].isin(valid_cmpd)).values
    result = run_pipeline(
        meta[ix],
        vals[ix],
        pos_sameby=['Metadata_JCP2022'],
        pos_diffby=[],
        neg_sameby=['Metadata_Plate'],
        neg_diffby=['Metadata_JCP2022'],
        null_size=10000,
        batch_size=20000,
        seed=0,
    )
    agg_result = aggregate(result, 'Metadata_JCP2022', threshold=0.05)
    agg_result = agg_result.query('Metadata_JCP2022 in @valid_cmpd')
    result.to_parquet(ap_path)
    agg_result.to_parquet(map_path)
