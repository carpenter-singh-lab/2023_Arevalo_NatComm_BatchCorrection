'''Correction methods'''
import logging

from harmonypy import run_harmony
import quality_control as qc

logger = logging.getLogger(__name__)


def harmony(dframe_path, batch_key, output_path):
    '''Harmony correction'''
    meta, feats, features = qc.io.split_parquet(dframe_path)
    # n_latent = min(feats.shape) - 1  # required for arpack
    # logger.info('Computing PCA...')
    # feats = sc.tl.pca(feats, n_comps=n_latent)  # Generates X_pca
    # logger.info('Computing PCA Done.')
    harmony_out = run_harmony(feats,
                              meta,
                              batch_key,
                              max_iter_harmony=20,
                              nclust=300)  # Number of compounds

    feats = harmony_out.Z_corr.T
    features = [f'harmony_{i}' for i in range(feats.shape[1])]
    qc.io.merge_parquet(meta, feats, features, output_path)
