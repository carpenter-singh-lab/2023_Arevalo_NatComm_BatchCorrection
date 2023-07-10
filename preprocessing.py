'''Functions to preprocess features'''
import logging
import anndata as ad
import numpy as np
import pandas as pd
from pycytominer.feature_select import feature_select

from cleaning import drop_na_inf, drop_outlier_feats, drop_outlier_samples
from loader import load_plates
from normalization import sphere
from utils import PathLocator, find_feat_cols, find_meta_cols, hashname

logger = logging.getLogger(__name__)


def filter_dmso(adata: ad.AnnData) -> ad.AnnData:
    '''Filter dmso'''
    dmso_idx = adata.obs.Metadata_JCP2022 == 'DMSO'
    adata = adata[~dmso_idx]
    return adata


def filter_values_norm(adata: ad.AnnData, locator: PathLocator):
    '''Filter values used for normalization'''
    meta = adata.obs
    field = locator.config['column_norm']
    values = locator.config['values_norm']
    norm_index = meta[field].isin(values)
    adata = adata[~norm_index]
    return adata


def filter_low_replicates(compound_data: ad.AnnData, min_replicates: int,
                          label_key: str):
    '''
    Filter compounds with min number of replicates
    '''
    counts = compound_data.obs[label_key].value_counts()
    logger.info(f'Before filter low replicate compounds: {counts.sum()}')
    counts = counts[(counts >= min_replicates)]
    logger.info(f'After filter low replicate compounds: {counts.sum()}')
    compound_data = compound_data[compound_data.obs[label_key].isin(
        counts.index)]
    return compound_data


def outlier_removal(dframe: pd.DataFrame):
    '''Remove outliers'''
    dframe, _ = drop_outlier_feats(dframe, threshold=1e2)
    dframe = drop_outlier_samples(dframe, threshold=1e2)
    # dframe = isolation_removal(dframe)
    return dframe


def feature_selection(dframe: pd.DataFrame):
    '''Run feature selection'''
    logger.info(f'{dframe.shape} before feature select')
    operations = [
        "variance_threshold", "correlation_threshold", "drop_na_columns",
        "blocklist"
    ]
    dframe = feature_select(dframe, operation=operations, image_features=True)
    logger.info(f'{dframe.shape} after feature select.')
    dframe = drop_na_inf(dframe, axis=1)
    return dframe


def to_anndata(dframe: pd.DataFrame):
    '''Convert pandas dataframe to AnnData object'''
    dframe.reset_index(drop=True, inplace=True)
    meta = dframe[find_meta_cols(dframe)]
    feat_cols = find_feat_cols(dframe)
    data = dframe[feat_cols].values.astype(np.float32)

    adata = ad.AnnData(data)
    adata.var_names = feat_cols
    adata.obs = meta
    return adata


def add_row_col(meta: pd.DataFrame):
    '''Add Metadata_Row and Metadata_Column to the DataFrame'''
    well_regex = r'^(?P<row>[a-zA-Z]{1,2})(?P<column>[0-9]{1,2})$'
    position = meta['Metadata_Well'].str.extract(well_regex)
    meta['Metadata_Row'] = position['row'].astype('category')
    meta['Metadata_Column'] = position['column'].astype('category')


def load_data(locator: PathLocator) -> ad.AnnData:
    '''Load data and saved processed data when needed'''
    config = locator.config
    hashid = hashname(config)
    logger.info(hashid)
    # Sources
    if locator.profiles_path.exists():
        logger.info('Loading preprocessed data...')
        adata = ad.read(locator.profiles_path)
        logger.info('Done.')
        add_row_col(adata.obs)
        return adata

    sources = config['sources']
    plate_types = config['plate_types']

    # Normalization per plate
    epsilon_mad = config['epsilon_mad']
    mad_norm = config['mad_norm']
    column_norm = config['column_norm']
    values_norm = config['values_norm']

    if locator.processed_plates_path.exists():
        dframe = pd.read_parquet(locator.processed_plates_path)
    else:
        logger.info('Loading plates...')
        dframe = load_plates(sources, plate_types, mad_norm, epsilon_mad,
                             column_norm, values_norm)
        logger.info('Loading plates done.')

        if config['outlier_removal']:
            dframe = outlier_removal(dframe)
        if config['nan_removal']:
            dframe = drop_na_inf(dframe)
        if config['feature_selection']:
            dframe = feature_selection(dframe)
        dframe = dframe.reset_index(drop=True).copy()
        dframe.to_parquet(locator.processed_plates_path)

    adata = to_anndata(dframe)
    adata.layers['raw'] = dframe[find_feat_cols(dframe)].values

    if config['sphering']:
        dframe, spherer = sphere(dframe, config['sphering_lambda'],
                                 config['sphering_mode'], column_norm,
                                 values_norm)
        adata.X = dframe[find_feat_cols(dframe)].values.astype('float32')
        np.savez_compressed(locator.spherer_path, spherer=spherer)
    adata.write(locator.profiles_path, compression='gzip')
    add_row_col(adata.obs)
    return adata
