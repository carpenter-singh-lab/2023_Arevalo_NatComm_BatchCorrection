from copairs.map import aggregate, run_pipeline
import pandas as pd

from quality_control.io import split_parquet


def _index(meta,
           plate_types,
           min_replicates,
           max_replicates,
           ignore_dmso=True):
    '''Select compounds to be used in mAP computation'''
    index = meta['Metadata_PlateType'].isin(plate_types)
    valid_cmpd = meta.loc[index, 'Metadata_JCP2022'].value_counts()
    if ignore_dmso:
        valid_cmpd.drop('DMSO', inplace=True)
    valid_cmpd = valid_cmpd[valid_cmpd.between(min_replicates, max_replicates)]
    valid_cmpd = valid_cmpd.index
    index = (index & meta['Metadata_JCP2022'].isin(valid_cmpd)).values
    return index


def average_precision(parquet_path,
                      ap_path,
                      plate_types,
                      min_replicates,
                      max_replicates,
                      ignore_dmso=True):
    meta, vals, _ = split_parquet(parquet_path)
    ix = _index(meta, plate_types, min_replicates, max_replicates, ignore_dmso)
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
    result.to_parquet(ap_path)


def mean_average_precision(ap_path, map_path, threshold=0.05):
    result = pd.read_parquet(ap_path)
    agg_result = aggregate(result, 'Metadata_JCP2022', threshold=threshold)
    agg_result.to_parquet(map_path)
