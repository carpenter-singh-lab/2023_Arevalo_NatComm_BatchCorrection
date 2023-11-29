'''
Functions to read results.
'''
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

from utils import METRIC_MAP, PathLocator, preprocessing_model_name

logger = logging.getLogger(__name__)


def load_scores(locator: PathLocator,
                redlist: list | None = None) -> pd.DataFrame:
    '''Load scores from a given locator in a pd.DataFrame object'''
    try:
        scores = pd.read_csv(locator.score_path, index_col='metric').score
        scores.dropna(inplace=True)
    except FileNotFoundError:
        logger.warning(f'single-cell scores not found for {locator.hashid}.{locator.model}')
        scores = pd.Series(name='score', dtype=np.float32)
        scores.index.name = 'metric'

    try:
        for key, val in load_map_compound(locator).items():
            scores[key] = val
    except FileNotFoundError:
        logger.warning(f'mAP scores not found for {locator.hashid}.{locator.model}')

    scores = scores.reset_index()
    scores['dimension'] = scores['metric'].map(METRIC_MAP)
    scores['batch_key'] = locator.config['batch_key']
    scores['label_key'] = locator.config['label_key']
    scores['model'] = locator.model
    scores['sphering'] = locator.config['sphering']
    scores['mad_norm'] = locator.config['mad_norm']
    scores['config_id'] = locator.hashid
    for key, val in locator.config.items():
        scores[key] = [val] * len(scores)

    if not locator.config['sphering_lambda']:
        # Lambda was automatically assigned. Try to recover from npz file.
        try:
            npz_obj = np.load(locator.spherer_path, allow_pickle=True)
            reg_val = npz_obj['spherer'].item().regularization
            scores['sphering_lambda'] = reg_val
        except FileNotFoundError:
            pass

    scores['model'] = scores['model'].apply(lambda x: preprocessing_model_name(x, locator.config))

    if redlist:
        scores = scores[~scores['metric'].isin(redlist)]

    return scores


def load_from_locators(locators: list[PathLocator],
                       redlist: list[str] | None = None,
                       verbose=0) -> pd.DataFrame:
    '''
    Load results from a list of locator objects.
    '''
    if not redlist:
        redlist = []
    redlist.append('null_th_repl')
    scores = []
    for locator in locators:
        try:
            score = load_scores(locator, redlist)
            score['min_replicates_to_eval'] = locator.min_replicates_to_eval
            scores.append(score)
        except FileNotFoundError as exc:
            if verbose:
                logger.warning(f'Results not ready: {exc.filename}')
    scores = pd.concat(scores)
    return scores


def load_same_conf(config_path: Path,
                   min_replicates_to_eval: int,
                   redlist: list[str] | None = None):
    '''
    Load results from multiple models under the same config.
    '''
    locators = []
    output_path = config_path.parent.parent
    items = config_path.parent.iterdir()
    models = [item.stem for item in items if item.is_dir()]
    for model in models:
        locator = PathLocator(config_path, model, min_replicates_to_eval,
                              output_path)
        locators.append(locator)
    return load_from_locators(locators, redlist)


def load_same_model(output_path: str, model: str, min_replicates_to_eval: int):
    '''
    Load results from multiple config for the same model.
    '''
    locators = []
    for path in find_configs(output_path):
        locator = PathLocator(path, model, min_replicates_to_eval, output_path)
        locators.append(locator)
    scores = load_from_locators(locators)
    return scores


def find_configs(output_path: str):
    '''Find all config.json files in a given folder'''
    config_paths = []
    for path in Path(output_path).iterdir():
        config_path = path / 'config.json'
        if config_path.is_file():
            config_paths.append(config_path)
    return config_paths


def filter_low_variance_metrics(scores: pd.DataFrame, tol: float):
    '''
    tolerance threshold to ignore metrics that does not change across methods.
    '''
    relevant = scores.groupby('metric').score.std()[lambda x: x > tol]
    scores = scores[scores.metric.isin(relevant.index)]
    return scores


def find_all_locators(output_path: str):
    '''Find all config experiments and create a PathLocator for each'''
    locators = []
    for config_path in find_configs(output_path):
        for child in config_path.parent.iterdir():
            if child.is_dir():
                model = child.stem
                for nrepl in child.iterdir():
                    if match := re.match(r'replicates_(\d+)', nrepl.stem):
                        min_repl = int(match.group(1))
                        locators.append(
                            PathLocator(config_path, model, min_repl,
                                        output_path))
    return locators


def load_all_results(output_path: str):
    '''Load all results from a given folder in a DataFrame indexed by 'sources'
    and plate_types'''
    scores = load_from_locators(find_all_locators(output_path))
    return scores
