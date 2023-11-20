import logging

import pandas as pd
from pycytominer import feature_select

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def select_features(dframe_path, feat_selected_path):
    '''Run feature selection'''
    dframe = pd.read_parquet(dframe_path)
    operations = [
        "variance_threshold",
        "correlation_threshold",
        # The two below were included in the first exploration. Ignoring them now
        # "drop_na_columns", "blocklist"
    ]
    for op in operations:
        before = dframe.shape[1]
        dframe = feature_select(dframe, operation=[op], image_features=True)
        after = dframe.shape[1]
        logger.info(f'{before} >> {after} features after {op}')
    dframe.reset_index(drop=True).to_parquet(feat_selected_path)
