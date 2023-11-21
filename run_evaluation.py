'''
Evaluate corrected representations
'''
import argparse
import logging
logging.basicConfig(format='%(levelname)s:%(asctime)s:%(name)s:%(message)s',
                    level=logging.INFO)
import anndata as ad
import pandas as pd
import scib
from copairs.replicating import CorrelationTestResult, correlation_test

import utils
from preprocessing import filter_low_replicates

logger = logging.getLogger(__name__)

METRICS_SET = {
    'ALL': [
        'NMI_cluster/label', 'ARI_cluster/label', 'ASW_label',
        'ASW_label/batch', 'PCR_batch', 'isolated_label_F1',
        'isolated_label_silhouette', 'graph_conn', 'kBET', 'iLISI', 'cLISI',
    ],
    'FAST': [
        'NMI_cluster/label', 'ARI_cluster/label', 'ASW_label',
        'ASW_label/batch', 'PCR_batch', 'graph_conn', 'iLISI', 'cLISI',
    ],
    'COLUMN_AS_BATCH': [
        'NMI_cluster/label', 'ARI_cluster/label', 'ASW_label',
        'ASW_label/batch', 'PCR_batch', 'isolated_label_silhouette',
        'graph_conn', 'iLISI', 'cLISI'
    ],
}


def percent_replicating(adata: ad.AnnData, use_rep: str,
                        null_size: int) -> CorrelationTestResult:
    '''Calculate percent replicating'''
    groupby = ['Metadata_JCP2022']
    diffby = ['Metadata_Plate', 'Metadata_Well']
    if adata.obs['Metadata_Source'].nunique() > 1:
        diffby = ['Metadata_Source'] + diffby
    return correlation_test(adata.obsm[use_rep], adata.obs, groupby, diffby,
                            null_size)


def sc_scores(adata: ad.AnnData, adata_int: ad.AnnData, use_rep: str,
              batch_key: str, label_key: str,
              metrics: list[str]) -> pd.Series:
    '''Report metrics from single-cell library'''

    isolated_labels_asw_='isolated_label_silhouette' in metrics
    silhouette_ = 'ASW_label' in metrics or 'ASW_label/batch' in metrics or isolated_labels_asw_
    params = dict(isolated_labels_asw_=isolated_labels_asw_,
                  silhouette_=silhouette_,
                  graph_conn_='graph_conn' in metrics,
                  pcr_='PCR_batch' in metrics,
                  isolated_labels_f1_='isolated_label_F1' in metrics,
                  nmi_='NMI_cluster/label' in metrics,
                  ari_='ARI_cluster/label' in metrics,
                  kBET_='kBET' in metrics,
                  ilisi_='iLISI' in metrics,
                  clisi_='cLISI' in metrics,
                  hvg_score_=False,
                  cell_cycle_=False,
                  trajectory_=False)
    scores = scib.metrics.metrics(adata,
                                  adata_int,
                                  batch_key,
                                  label_key,
                                  embed=use_rep,
                                  **params)
    return scores[0]


def evaluate_sc(locator: utils.PathLocator, compound_data: ad.AnnData,
                corrected_data: ad.AnnData, metrics: list[str]):
    '''
    Compute and save scores for corrected data under this experiment
    '''
    model = locator.model
    use_rep = f'X_{model}'

    scores = sc_scores(compound_data, corrected_data, use_rep,
                       locator.config['batch_key'],
                       locator.config['label_key'], metrics)

    if locator.score_path.exists():
        logger.info(f'Single cell scores for {model} found.')
        prev_scores = pd.read_csv(locator.score_path)
        prev_scores = prev_scores.set_index('metric').score
        scores.update(prev_scores)

    scores = scores.reset_index()
    scores.columns = ['metric', 'score']
    scores.to_csv(locator.score_path, index=False)


def evaluate_percent_repl(locator: utils.PathLocator,
                          corrected_data: ad.AnnData):
    '''compute score based on correlation replicability'''
    model = locator.model
    use_rep = f'X_{model}'
    if locator.corr_dist_path.exists() and locator.null_dist_path.exists():
        logger.info(f'Percent replicating for {model} found.')
        return

    logger.info('Percent replicating...')
    corr_result = percent_replicating(corrected_data, use_rep, null_size=1000)
    corr_result.corr_df.to_csv(locator.corr_dist_path)
    corr_result.null_dist.to_csv(locator.null_dist_path, index=False)
    logger.info('Percent replicating done.')


def workflow(locator: utils.PathLocator, metrics: list[str]):
    '''Main workflow'''
    # try to load precomputed metrics
    try:
        precomputed = pd.read_csv(locator.score_path).dropna().metric.tolist()
    except FileNotFoundError:
        precomputed = []
    if locator.corr_dist_path.exists() and locator.null_dist_path.exists():
        precomputed.append('percent_repl')

    metrics = [metric for metric in metrics if metric not in precomputed]
    if not metrics:
        return

    # Load corrected data.
    logger.info(f'Computing metrics: {metrics}.')
    logger.info(f'Loading corrected data for {locator.model}...')
    corrected_data = ad.read(locator.corrected_path)
    logger.info('Load corrected data done.')

    # Filter low replicates before computing metrics
    corrected_data = filter_low_replicates(corrected_data,
                                           locator.min_replicates_to_eval,
                                           locator.config['label_key'])

    # Use data from raw layer to compute PCR_batch.
    corrected_data.X = corrected_data.layers['raw']

    # Compute metrics
    evaluate_sc(locator, corrected_data, corrected_data, metrics)
    if 'percent_repl' in metrics:
        evaluate_percent_repl(locator, corrected_data)


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
        description=('Compute evaluation metrics following the config file'))
    parser.add_argument('config_path', type=str)
    parser.add_argument(
        'model',
        help='Model to evaluate. "all" to evaluate all the correction methods',
        choices=models + ['all'])
    parser.add_argument('min_replicates_to_eval', type=int)
    parser.add_argument('output_path', type=str)
    parser.add_argument(
        'metrics',
        nargs='*',
        default=['ALL'],
        help=('Metrics to compute. Default ALL. Valid values are '
              'ALL: All metrics, or '
              'FAST: Only metrics that are not computationally expensive, or '
              'a custom list of metrics'))

    args = parser.parse_args()
    if 'ALL' in args.metrics:
        metrics = METRICS_SET['ALL']
    elif 'FAST' in args.metrics:
        metrics = METRICS_SET['FAST']
    elif 'COLUMN_AS_BATCH' in args.metrics:
        metrics = METRICS_SET['COLUMN_AS_BATCH']
    else:
        metrics = args.metrics

    if args.model == 'all':
        for model in models:
            locator = utils.PathLocator(args.config_path, model,
                                        args.min_replicates_to_eval,
                                        args.output_path)
            workflow(locator, metrics)

    else:
        locator = utils.PathLocator(args.config_path, args.model,
                                    args.min_replicates_to_eval,
                                    args.output_path)
        workflow(locator, metrics)


if __name__ == "__main__":
    main()
