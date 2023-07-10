'''
Evaluate compound replicability using mAP
'''
import argparse
import logging

import anndata as ad
from copairs.map import run_pipeline

from utils import PathLocator

logger = logging.getLogger(__name__)


def evaluate_map(locator: PathLocator, corrected_data: ad.AnnData):
    '''compute score based on mean average precision'''
    null_size, batch_size = 10000, 100000
    meta = corrected_data.obs
    feats = corrected_data.obsm[f'X_{locator.model}']

    pos_sameby = neg_diffby = locator.config['label_key']
    pos_diffby = locator.config['batch_key']
    neg_sameby = 'Metadata_Plate'
    result = run_pipeline(
        meta,
        feats,
        pos_sameby,
        pos_diffby,
        neg_sameby,
        neg_diffby,
        null_size,
        batch_size=batch_size,
    )

    result.to_csv(locator.average_precision_path, index=False)

    return result


def workflow(locator: PathLocator):
    # Load corrected data.
    if locator.average_precision_path.is_file():
        return
    logger.info(f'Loading corrected data for {locator.model}...')
    corrected_data = ad.read(locator.corrected_path)
    logger.info('Load corrected data done.')
    evaluate_map(locator, corrected_data)


def main():
    '''Parse input params'''
    models = [
        'sphering',
        'combat',
        'harmony',
        'scanorama',
        'mnn',
        'desc',
        'scvi',
    ]

    parser = argparse.ArgumentParser(
        description=('Compute map for replicates'))
    parser.add_argument('config_path', type=str)
    parser.add_argument(
        'model',
        help='Model to evaluate. "all" to evaluate all the correction methods',
        choices=models + ['all'])
    parser.add_argument('output_path', type=str)

    args = parser.parse_args()

    if args.model == 'all':
        for model in models:
            locator = PathLocator(args.config_path, model, 5, args.output_path)
            workflow(locator)

    else:
        locator = PathLocator(args.config_path, args.model, 5,
                              args.output_path)
        workflow(locator)


if __name__ == "__main__":
    main()
