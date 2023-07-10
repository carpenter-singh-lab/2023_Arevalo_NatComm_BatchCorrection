'''
Functions to browse the folders and load information
'''
from functools import partial
import logging

import numpy as np
import pandas as pd
from pycytominer.operations.transform import RobustMAD
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map

from normalization import apply_norm
from utils import FEATURE_SET

logger = logging.getLogger(__name__)

MAPPER = {
    'JCP2022_085227': 'Aloxistatin',
    'JCP2022_037716': 'AMG900',
    'JCP2022_025848': 'Dexamethasone',
    'JCP2022_046054': 'FK-866',
    'JCP2022_035095': 'LY2109761',
    'JCP2022_064022': 'NVS-PAK1-1',
    'JCP2022_050797': 'Quinidine',
    'JCP2022_012818': 'TC-S-7004',
    'JCP2022_033924': 'DMSO',
    'JCP2022_999999': 'UNTREATED',
    'JCP2022_UNKNOWN': 'UNKNOWN',
    'JCP2022_900001': 'BAD CONSTRUCT'
}

MICRO_CONFIG = pd.read_csv('https://raw.githubusercontent.com/jump-cellpainting/datasets/181fa0dc96b0d68511b437cf75a712ec782576aa/metadata/microscope_config.csv')
MICRO_CONFIG['Metadata_Source'] = 'source_' + MICRO_CONFIG['Metadata_Source'].astype(str)
MICRO_CONFIG = MICRO_CONFIG.set_index('Metadata_Source')['Metadata_Microscope_Name']

def get_source_4_plate_redlist(plate_types: list[str]):
    '''Get set of plate_id's  that should be not considered in the analysis'''
    # https://github.com/jump-cellpainting/jump-orf-analysis/issues/1#issuecomment-921888625
    # Low concentration plates
    redlist = set(['BR00127147', 'BR00127148', 'BR00127145', 'BR00127146'])
    # https://github.com/jump-cellpainting/aws/issues/70#issuecomment-1182444836
    redlist.add('BR00123528A')

    metadata = pd.read_csv('inputs/experiment-metadata.tsv', sep='\t')
    if 'ORF' in plate_types:
        # filter ORF plates.
        query = 'Batch=="Batch12"'
        bad_plates = set(metadata.query(query).Assay_Plate_Barcode)
        redlist |= bad_plates

    if 'TARGET2' in plate_types:
        # filter TARGET2 plates
        query = 'Anomaly!="none"'
        bad_plates = set(metadata.query(query).Assay_Plate_Barcode)
        redlist |= bad_plates
    return redlist


SOURCE3_BATCH_REDLIST = {
    'CP_32_all_Phenix1', 'CP_33_all_Phenix1', 'CP_34_mix_Phenix1',
    'CP_35_all_Phenix1', 'CP_36_all_Phenix1', 'CP59', 'CP60'
}


def build_path(row: pd.Series) -> str:
    '''Create the path to the parquet file'''
    template = ('./inputs/{Metadata_Source}/workspace/profiles/'
                '{Metadata_Batch}/{Metadata_Plate}/{Metadata_Plate}.parquet')
    return template.format(**row.to_dict())


def get_plate_metadata(sources: list[str],
                       plate_types: list[str]) -> pd.DataFrame:
    '''Create filtered metadata DataFrame'''
    plate_metadata = pd.read_csv('./inputs/metadata/plate.csv.gz')
    # Filter plates from source_4
    if 'source_4' in sources:
        redlist = get_source_4_plate_redlist(plate_types)
        plate_metadata = plate_metadata[~plate_metadata['Metadata_Plate'].
                                        isin(redlist)]

    # Filter plates from source_3 batches without DMSO
    plate_metadata = plate_metadata[
        (~plate_metadata['Metadata_Batch'].isin(SOURCE3_BATCH_REDLIST))
        | (plate_metadata['Metadata_PlateType'] == 'TARGET2')]

    plate_metadata = plate_metadata[plate_metadata['Metadata_Source'].isin(
        sources)]
    plate_metadata = plate_metadata[plate_metadata['Metadata_PlateType'].isin(
        plate_types)]
    return plate_metadata


def get_well_metadata(plate_types: list[str]):
    '''Load well metadata'''
    well_metadata = pd.read_csv('./inputs/metadata/well.csv.gz')
    if 'ORF' in plate_types:
        orf_metadata = pd.read_csv('./inputs/metadata/orf.csv.gz')
        well_metadata = well_metadata.merge(orf_metadata, how='inner')
    # Use readable names for controls and non-treatment codes
    well_metadata['Metadata_JCP2022'] = well_metadata[
        'Metadata_JCP2022'].apply(lambda x: MAPPER.get(x, x))
    # Filter out wells
    well_metadata = well_metadata[~well_metadata['Metadata_JCP2022'].isin(
        ['UNTREATED', 'UNKNOWN', 'BAD CONSTRUCT'])]

    return well_metadata


def load_plate(path: str, meta: pd.DataFrame, mad_norm: bool,
               epsilon_mad: float, column_norm: str,
               values_norm: list[str]) -> pd.DataFrame:
    '''Load a plate and optionally apply robustmad norm'''
    profile = pd.read_parquet(path)
    for column in FEATURE_SET:
        profile[column] = profile[column].astype(np.float32)
    profile['Metadata_Plate'] = profile['Metadata_Plate'].astype(str)
    plate = profile.merge(
        meta, on=['Metadata_Source', 'Metadata_Plate', 'Metadata_Well'])
    if plate.shape[0] == 0:
        logger.warning(f'{path} plate gives empty result')
        return plate

    if mad_norm:
        normalizer = RobustMAD(epsilon_mad)
        plate = apply_norm(normalizer, plate, column_norm, values_norm)
    return plate


def filter_compounds_no_dmso(meta: pd.DataFrame):
    '''Filter out compound plates without any DMSO'''
    compound = meta.query('Metadata_PlateType=="COMPOUND"')
    plate_ids = compound['Metadata_Plate'].drop_duplicates()
    plate_with_dmso = compound.query(
        'Metadata_JCP2022=="DMSO"').Metadata_Plate.drop_duplicates()
    ignored_plates = plate_ids[~plate_ids.isin(plate_with_dmso)]
    if len(ignored_plates) > 0:
        logger.warning(
            f'Ignoring {ignored_plates.values} compound plates, do not have DMSO wells'
        )
        meta = meta[~meta['Metadata_Plate'].isin(ignored_plates)]
    return meta

def load_metadata(sources: list[str], plate_types: list[str]):
    '''Load metadata only'''
    plate = get_plate_metadata(sources, plate_types)
    well = get_well_metadata(plate_types)
    meta = well.merge(plate, on=['Metadata_Source', 'Metadata_Plate'])
    return meta

def load_plates(sources: list[str],
                plate_types: list[str],
                mad_norm: bool,
                epsilon_mad: float,
                column_norm: str,
                values_norm: list[str],
                parallel: bool = True) -> pd.DataFrame:
    '''Load plates'''
    plate_metadata = get_plate_metadata(sources, plate_types)
    well_metadata = get_well_metadata(plate_types)
    paths = plate_metadata.apply(build_path, axis=1).values

    meta = well_metadata.merge(plate_metadata,
                               on=['Metadata_Source', 'Metadata_Plate'])
    meta = filter_compounds_no_dmso(meta)
    par_func = partial(load_plate,
                       meta=meta,
                       mad_norm=mad_norm,
                       epsilon_mad=epsilon_mad,
                       column_norm=column_norm,
                       values_norm=values_norm)
    if parallel:
        plates = process_map(par_func, paths, chunksize=1)
    else:
        plates = map(par_func, tqdm(paths))
    dframe = pd.concat(plates)
    for col in dframe:
        if col not in FEATURE_SET:
            dframe[col] = dframe[col].astype('category')
    return dframe
