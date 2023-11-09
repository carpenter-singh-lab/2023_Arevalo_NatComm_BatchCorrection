import json

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from tqdm.contrib.concurrent import thread_map

from loader import load_metadata, build_path
from utils import find_feat_cols, find_meta_cols


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
    paths = (meta[['Metadata_Source', 'Metadata_Batch', 'Metadata_Plate'
                   ]].drop_duplicates().apply(build_path, axis=1)).values

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


def write_data(config_path, output_file):
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
