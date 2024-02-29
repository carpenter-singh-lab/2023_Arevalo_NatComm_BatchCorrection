import logging
import os
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import silhouette_score

from preprocessing.io import split_parquet, to_anndata

with warnings.catch_warnings():
    warnings.filterwarnings('ignore',
                            category=ResourceWarning,
                            message='Implicitly cleaning up')
    from scib import metrics
    from scib.metrics.pcr import pc_regression

logger = logging.getLogger(__name__)
CLUSTER_KEY = 'Metadata_Cluster'


def filter_dmso_anndata(parquet_path):
    adata = to_anndata(parquet_path)
    non_dmso_ix = adata.obs['Metadata_JCP2022'] != 'DMSO'
    return adata[non_dmso_ix].copy()


def filter_dmso(parquet_path):
    meta, feats, features = split_parquet(parquet_path)
    non_dmso_ix = meta['Metadata_JCP2022'] != 'DMSO'
    meta = meta[non_dmso_ix].reset_index(drop=True).copy()
    feats = feats[non_dmso_ix]
    return meta, feats, features


def cluster(parquet_path, adata_path):
    adata = filter_dmso_anndata(parquet_path)
    logger.info('compute neighbors')
    sc.pp.neighbors(adata, use_rep='X', n_neighbors=25, metric='cosine')
    logger.info('run clustering')
    sc.tl.leiden(adata, key_added=CLUSTER_KEY)
    adata.write_h5ad(adata_path, compression='gzip')




def nmi(adata_path, label_key, nmi_path):
    adata = ad.read_h5ad(adata_path)
    nmi = metrics.nmi(adata, label_key, CLUSTER_KEY)
    np.array(nmi).tofile(nmi_path)


def ari(adata_path, label_key, ari_path):
    adata = ad.read_h5ad(adata_path)
    ari = metrics.ari(adata, label_key, CLUSTER_KEY)
    np.array(ari).tofile(ari_path)


def asw(parquet_path, label_key, asw_path):
    meta, feats, _ = filter_dmso(parquet_path)
    asw = silhouette_score(feats, meta[label_key], metric='cosine')
    asw = (asw + 1) / 2
    np.array(asw).tofile(asw_path)


def silhouette_batch(parquet_path, label_key, batch_key, asw_batch_path):
    adata = filter_dmso_anndata(parquet_path)
    adata.obsm['X_embd'] = adata.X
    asw_batch = metrics.silhouette_batch(
        adata,
        batch_key,
        label_key,
        'X_embd',
        metric="cosine",
        verbose=False,
    )
    np.array(asw_batch).tofile(asw_batch_path)


def pcr_batch(pre_parquet_path, post_parquet_path, batch_key, pcr_batch_path):
    adata = filter_dmso_anndata(pre_parquet_path)
    adata_int = filter_dmso_anndata(post_parquet_path)
    pcr_score = metrics.pcr_comparison(adata,
                                       adata_int,
                                       embed=None,
                                       covariate=batch_key,
                                       verbose=False)
    np.array(pcr_score).tofile(pcr_batch_path)


def pcr(parquet_path, batch_key, pcr_path):
    meta, vals, _ = filter_dmso(parquet_path)
    # 1 - pcr to make it 0 worst 1 best
    pcr_score = 1 - pc_regression(vals, meta[batch_key].values)
    np.array(pcr_score).tofile(pcr_path)


def isolated_labels_asw(adata_path, label_key, batch_key, il_asw_path):
    adata = ad.read_h5ad(adata_path)
    adata.obsm['X_embd'] = adata.X
    il_score_asw = metrics.isolated_labels(
        adata,
        label_key=label_key,
        batch_key=batch_key,
        embed='X_embd',
        cluster=False,
        # max number of batches a label should be present to be considered isolated
        iso_threshold=None,
        verbose=False,
    )
    np.array(il_score_asw).tofile(il_asw_path)


def isolated_labels_f1(adata_path, label_key, batch_key, il_f1_path):
    adata = ad.read_h5ad(adata_path)
    adata.obsm['X_embd'] = adata.X
    il_score_f1 = metrics.isolated_labels(
        adata,
        label_key=label_key,
        batch_key=batch_key,
        embed='X_embd',
        cluster=True,
        # max number of batches a label should be present to be considered isolated
        iso_threshold=None,
        verbose=False,
    )
    np.array(il_score_f1).tofile(il_f1_path)


def graph_connectivity(adata_path, label_key, graph_conn_path):
    adata = ad.read_h5ad(adata_path)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',
                                category=FutureWarning,
                                message='pandas.value_counts is deprecated')
        graph_conn_score = metrics.graph_connectivity(adata, label_key=label_key)
        np.array(graph_conn_score).tofile(graph_conn_path)


def kbet(adata_path, label_key, batch_key, kbet_path):
    adata = ad.read_h5ad(adata_path)
    M = metrics.kbet.diffusion_conn(adata, min_k=15, copy=False)
    adata.obsp["connectivities"] = M
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', module='rpy2.robjects')
        warnings.filterwarnings('ignore', category=FutureWarning)
        kbet_score = metrics.kBET(
            adata,
            batch_key=batch_key,
            label_key=label_key,
            type_='knn',
            embed=None,
            scaled=True,
            verbose=False,
        )
        np.array(kbet_score).tofile(kbet_path)


def lisi_label(adata_path, label_key, lisi_label_path):
    adata = ad.read_h5ad(adata_path)
    clisi = metrics.clisi_graph(
        adata,
        label_key=label_key,
        type_='knn',
        subsample=100,  # Use all data
        scale=True,
        n_cores=8,
        verbose=False,
    )
    np.array(clisi).tofile(lisi_label_path)


def lisi_batch(adata_path, batch_key, lisi_batch_path):
    adata = ad.read_h5ad(adata_path)
    ilisi = metrics.ilisi_graph(
        adata,
        batch_key=batch_key,
        type_='knn',
        subsample=100,  # Use all data
        scale=True,
        n_cores=8,
        verbose=False,
    )
    np.array(ilisi).tofile(lisi_batch_path)


def concat(*metric_paths, output_path):
    '''Concatenate scib metrics in a single file'''
    # Extract metric names from path
    start = len(os.path.commonprefix(metric_paths))
    end = -len('.bin')
    metrics = map(lambda x: x[start:end], metric_paths)

    # Concat metric values
    scores = np.fromiter(map(np.fromfile, metric_paths), dtype=np.float32)
    dframe = pd.DataFrame({'metric': metrics, 'score': scores})
    dframe.to_parquet(output_path)
