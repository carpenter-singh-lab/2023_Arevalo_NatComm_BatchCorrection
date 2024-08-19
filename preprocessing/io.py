import numpy as np
import pandas as pd
# import fastparquet
import pyarrow.parquet as pq
from tqdm.contrib.concurrent import thread_map
import anndata as ad

from .metadata import build_path, load_metadata, MICRO_CONFIG, find_feat_cols, find_meta_cols


def to_anndata(parquet_path):
    meta, feats, features = split_parquet(parquet_path)
    meta.index = meta.index.astype(str)
    adata = ad.AnnData(feats, meta)
    adata.var_names = features
    return adata


def split_parquet(dframe_path,
                  features=None) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    dframe = pd.read_parquet(dframe_path)
    if features is None:
        features = find_feat_cols(dframe)
    vals = np.empty((len(dframe), len(features)), dtype=np.float32)
    for i, c in enumerate(features):
        vals[:, i] = dframe[c]
    meta = dframe[find_meta_cols(dframe)].copy()
    return meta, vals, features


def merge_parquet(meta, vals, features, output_path) -> None:
    '''Save the data in a parquet file resetting the index'''
    dframe = pd.DataFrame(vals, columns=features)
    for c in meta:
        dframe[c] = meta[c].reset_index(drop=True)
    dframe.to_parquet(output_path)

# def get_num_rows(path) -> int:
#     '''Count the number of rows in a parquet file using fastparquet'''
#     pf = fastparquet.ParquetFile(path)
#     return pf.count()

def get_num_rows(path) -> int:
    '''Count the number of rows in a parquet file'''
    with pq.ParquetFile(path) as file:
        return file.metadata.num_rows


def prealloc_params(sources, plate_types):
    '''
    Get a list of paths to the parquet files and the corresponding slices
    for further concatenation
    '''
    meta = load_metadata(sources, plate_types)
    paths = (meta[['Metadata_Source', 'Metadata_Batch',
                   'Metadata_Plate']].drop_duplicates().apply(build_path,
                                                              axis=1)).values

    counts = thread_map(get_num_rows, paths, leave=False, desc='counts')
    slices = np.zeros((len(paths), 2), dtype=int)
    slices[:, 1] = np.cumsum(counts)
    slices[1:, 0] = slices[:-1, 1]
    return paths, slices


def load_data(sources, plate_types):
    '''Load all plates given the params'''
    paths, slices = prealloc_params(sources, plate_types)
    total = slices[-1, 1]

    # f = fastparquet.ParquetFile(paths[0])
    # meta_cols = find_meta_cols(f.columns)
    # feat_cols = find_feat_cols(f.columns)

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
    for col in meta_cols:
        dframe[col] = meta[col]
    return dframe


def add_pert_type(meta: pd.DataFrame, col: str = 'Metadata_PertType'):
    meta[col] = 'trt'
    meta.loc[~meta['Metadata_JCP2022'].str.startswith('JCP'), col] = 'poscon'
    meta.loc[meta['Metadata_JCP2022'] == 'DMSO', col] = 'negcon'
    meta[col] = meta[col].astype('category')


def add_row_col(meta: pd.DataFrame):
    '''Add Metadata_Row and Metadata_Column to the DataFrame'''
    well_regex = r'^(?P<row>[a-zA-Z]{1,2})(?P<column>[0-9]{1,2})$'
    position = meta['Metadata_Well'].str.extract(well_regex)
    meta['Metadata_Row'] = position['row'].astype('category')
    meta['Metadata_Column'] = position['column'].astype('category')

def add_microscopy_info(meta: pd.DataFrame):
    configs = meta['Metadata_Source'].map(MICRO_CONFIG).astype('category')
    meta['Metadata_Microscope'] = configs

def write_parquet(sources, plate_types, output_file):
    '''Write the parquet dataset given the params'''
    dframe = load_data(sources, plate_types)
    # Efficient merge
    meta = load_metadata(sources, plate_types)
    add_pert_type(meta)
    add_row_col(meta)
    add_microscopy_info(meta)
    foreign_key = ['Metadata_Source', 'Metadata_Plate', 'Metadata_Well']
    meta = dframe[foreign_key].merge(meta, on=foreign_key, how='left')
    for c in meta:
        dframe[c] = meta[c].astype('category')
    # Dropping samples with no metadata
    dframe.dropna(subset=['Metadata_JCP2022'], inplace=True)
    dframe.reset_index(drop=True, inplace=True)
    dframe.to_parquet(output_file)
