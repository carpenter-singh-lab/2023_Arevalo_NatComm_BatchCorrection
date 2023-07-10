'''
Run correction algorithms
'''
import argparse
import logging
from typing import Literal

import anndata as ad
import pandas as pd
import scanpy as sc
from scvi.model.utils import mde

import correct
from preprocessing import filter_values_norm, load_data
from utils import PathLocator

logging.basicConfig(format='%(levelname)s:%(asctime)s:%(name)s:%(message)s',
                    level=logging.INFO)
logger = logging.getLogger(__name__)


def obsm_to_df(adata: ad.AnnData, obsm_key: str,
               columns: list[str] | None) -> pd.DataFrame:
    '''Convert AnnData object to DataFrame using obs and obsm properties'''
    meta = adata.obs.reset_index(drop=True)
    feats = adata.obsm[obsm_key]
    n_feats = feats.shape[1]
    if not columns:
        columns = [f'{obsm_key}_{i:04d}' for i in range(n_feats)]
    data = pd.DataFrame(feats, columns=columns)
    return pd.concat([meta, data], axis=1)


def embd_2d(locator: PathLocator, corrected_data: ad.AnnData,
            vis: Literal['mde', 'pca', 'umap']):
    vmap = {
        'umap': locator.umap_path,
        'mde': locator.mde_path,
        'pca': locator.pca_path
    }
    path = vmap[vis]
    if path.exists():
        logger.info(f'{vis} embeddings found.')
        return

    use_rep = f'X_{locator.model}'
    if vis == 'mde':
        embd = corrected_data
        embd.obsm['X_mde'] = mde(corrected_data.obsm[use_rep])
    elif vis == 'pca':
        embd = ad.AnnData(corrected_data.obsm[use_rep], corrected_data.obs)
        sc.tl.pca(embd, n_comps=2)  # Generates X_pca
    elif vis == 'umap':
        use_rep = f'X_{locator.model}'
        embd = ad.AnnData(corrected_data.obsm[use_rep], corrected_data.obs)
        sc.pp.neighbors(embd, use_rep='X')
        sc.tl.umap(embd)  # Generates X_umap
    else:
        raise ValueError(f'{vis}: Invalid visualization method')
    dframe = obsm_to_df(embd, f'X_{vis}', ['x', 'y'])
    dframe.to_csv(path, index=False)


def run_correction(locator: PathLocator, compound_data: ad.AnnData):
    '''Run correction model in the compound data'''
    model = locator.model

    if locator.corrected_path.exists():
        logger.info(f'Loading corrected data for {model}...')
        compound_data = ad.read(locator.corrected_path)
        logger.info('Load corrected data done.')
        return compound_data

    logger.info(f'Running {model} with config {locator.hashid}')
    use_rep = f'X_{model}'
    correction_map = correct.get_method_map(locator.config['batch_key'],
                                            locator.config['label_key'])
    correct_fn = correction_map[model]
    correct_fn(compound_data, corrected_embed=use_rep)
    compound_data.X = compound_data.X.astype('float32')
    compound_data.write(locator.corrected_path, compression='gzip')
    return compound_data


def workflow(locator: PathLocator, vis: Literal['pca', 'umap', 'mde']):
    '''Main workflow'''
    vmap = {
        'umap': locator.umap_path,
        'mde': locator.mde_path,
        'pca': locator.pca_path
    }
    embd_path = vmap[vis]
    if embd_path.exists() and locator.corrected_path.exists():
        return
    # Load data
    full_data = load_data(locator)
    # Remove dmso
    compound_data = filter_values_norm(full_data, locator).copy()
    del full_data
    # Sort samples per batch_key. Required for scanorama
    batch_col = compound_data.obs[locator.config['batch_key']]
    idx = batch_col.argsort().values
    compound_data = compound_data[idx, :]
    # Correct data
    corrected_data = run_correction(locator, compound_data.copy())
    # Compute 2D embeddings
    embd_2d(locator, corrected_data, vis)


def main():
    '''Parse input params'''
    parser = argparse.ArgumentParser(description=(
        'Run batch correction algorithm following the config file'), )
    models = [
        'sphering',
        'combat',
        'harmony',
        'scanorama',
        'mnn',
        'desc',
        'scvi',
    ]
    parser.add_argument('config_path',
                        type=str,
                        help='json file with params. See config.json')
    parser.add_argument(
        'model',
        help='Model to correct data. "all" will run all the correction methods',
        choices=models + ['all'])
    parser.add_argument('output_path', type=str)

    parser.add_argument(
        '--vis',
        default='mde',
        help='Visualization embedding method. "pca", "umap" or "mde"',
        choices=['pca', 'umap', 'mde'])

    args = parser.parse_args()
    if args.model == 'all':
        for model in models:
            locator = PathLocator(
                args.config_path,
                model,
                5,  # TODO: Decouple this param from PathLocator
                args.output_path)
            workflow(locator, args.vis)
    else:
        locator = PathLocator(args.config_path, args.model, 5,
                              args.output_path)
        workflow(locator, args.vis)


if __name__ == "__main__":
    main()
