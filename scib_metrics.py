import os
import logging
import warnings

import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
from scib import metrics
from sklearn.metrics import silhouette_score

from quality_control.io import split_parquet, to_anndata

logger = logging.getLogger(__name__)

CLUSTER_KEY = 'Metadata_Cluster'


def cluster(parquet_path, label_key, adata_path):
    adata = to_anndata(parquet_path)
    logger.info('compute neighbors')
    sc.pp.neighbors(adata, use_rep='X', n_neighbors=15, metric='cosine')
    # Get cluster and neighbors
    logger.info('run clustering')
    metrics.cluster_optimal_resolution(adata,
                                       label_key=label_key,
                                       cluster_key=CLUSTER_KEY,
                                       metric=metrics.nmi,
                                       verbose=False)
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
    meta, feats, _ = split_parquet(parquet_path)
    asw = silhouette_score(feats, meta[label_key], metric='cosine')
    asw = (asw + 1) / 2
    np.array(asw).tofile(asw_path)


def silhouette_batch(parquet_path, label_key, batch_key, asw_batch_path):
    adata = to_anndata(parquet_path)
    adata.obsm['X_sphe'] = adata.X
    asw_batch = metrics.silhouette_batch(
        adata,
        batch_key,
        label_key,
        'X_sphe',
        metric="cosine",
        verbose=False,
    )
    np.array(asw_batch).tofile(asw_batch_path)


def pcr_batch(pre_parquet_path, post_parquet_path, batch_key, pcr_batch_path):
    adata = to_anndata(pre_parquet_path)
    adata_int = to_anndata(post_parquet_path)
    pcr_score = metrics.pcr_comparison(adata,
                                       adata_int,
                                       embed=None,
                                       covariate=batch_key,
                                       verbose=False)
    np.array(pcr_score).tofile(pcr_batch_path)


def isolated_labels_asw(adata_path, label_key, batch_key, il_asw_path):
    adata = ad.read_h5ad(adata_path)
    adata.obsm['X_sphe'] = adata.X
    il_score_asw = metrics.isolated_labels(
        adata,
        label_key=label_key,
        batch_key=batch_key,
        embed='X_sphe',
        cluster=False,
        # max number of batches a label should be present to be considered isolated
        iso_threshold=None,
        verbose=False,
    )
    np.array(il_score_asw).tofile(il_asw_path)


def isolated_labels_f1(adata_path, label_key, batch_key, il_f1_path):
    adata = ad.read_h5ad(adata_path)
    adata.obsm['X_sphe'] = adata.X
    il_score_f1 = metrics.isolated_labels(
        adata,
        label_key=label_key,
        batch_key=batch_key,
        embed='X_sphe',
        cluster=True,
        # max number of batches a label should be present to be considered isolated
        iso_threshold=None,
        verbose=False,
    )
    np.array(il_score_f1).tofile(il_f1_path)


def graph_connectivity(adata_path, label_key, graph_conn_path):
    adata = ad.read_h5ad(adata_path)
    graph_conn_score = metrics.graph_connectivity(adata, label_key=label_key)
    np.array(graph_conn_score).tofile(graph_conn_path)


def kbet(adata_path, label_key, batch_key, kbet_path):
    adata = ad.read_h5ad(adata_path)
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
