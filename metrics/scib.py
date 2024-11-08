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
    warnings.filterwarnings(
        "ignore",
        category=ResourceWarning,
        message="Implicitly cleaning up"
    )
    import scib_metrics as metrics
    from scib_metrics.utils import principal_component_regression as pc_regression

logger = logging.getLogger(__name__)
CLUSTER_KEY = "Metadata_Cluster"


def _add_moa_info(meta) -> pd.DataFrame:
    # from https://github.com/theislab/jump-integrate-reproducibility/blob/main/experiments/get_moas_from_opentarget/get_moas_from_opentarget.ipynb
    moa_info = pd.read_parquet("inputs/metadata/opentargets_moa_target2_eval.parquet")

    return meta.merge(moa_info, on="Metadata_InChIKey", how="left")

def _subset_obs_based_on_eval_label_presence(meta, feats, eval_key) -> pd.DataFrame:
    # using mask because pd.DataFrames and np.arrays index differently
    valid_mask = meta[eval_key].isna().values
    feats_subset = feats[valid_mask]
    meta_subset = meta.loc[valid_mask].copy()
    
    return feats_subset, meta_subset


def filter_dmso_anndata(parquet_path):
    adata = to_anndata(parquet_path)
    non_dmso_ix = adata.obs["Metadata_JCP2022"] != "DMSO"
    return adata[non_dmso_ix].copy()


def filter_dmso(parquet_path):
    meta, feats, features = split_parquet(parquet_path)
    non_dmso_ix = meta["Metadata_JCP2022"] != "DMSO"
    meta = meta[non_dmso_ix].reset_index(drop=True).copy()
    feats = feats[non_dmso_ix]
    return meta, feats, features


def cluster(parquet_path, adata_path):
    adata = filter_dmso_anndata(parquet_path)
    logger.info("compute neighbors")
    sc.pp.neighbors(adata, use_rep="X", n_neighbors=25, metric="cosine")
    logger.info("run clustering")
    sc.tl.leiden(
        adata,
        key_added=CLUSTER_KEY, 
        # Next 3 parameters to silence FutureWarning
        flavor="igraph",
        n_iterations=2, 
        directed=False
    ) 
    adata.write_h5ad(adata_path, compression="gzip")


def nmi(adata_path, label_key, nmi_path):
    adata = ad.read_h5ad(adata_path)
    nmi = metrics.nmi(adata, label_key, CLUSTER_KEY)

    np.save(nmi_path, nmi)


def ari(adata_path, label_key, ari_path):
    adata = ad.read_h5ad(adata_path)
    ari = metrics.ari(adata, label_key, CLUSTER_KEY)

    np.save(ari_path, ari)


def asw(parquet_path, eval_key, asw_path):

    meta, feats, _ = filter_dmso(parquet_path)
    meta = _add_moa_info(meta)

    if eval_key not in meta.columns:
        raise ValueError(f"Eval key '{eval_key}' not in metadata")

    meta = meta.reset_index(drop=True)
    feats_df = pd.DataFrame(feats)
    feats_df = feats_df.reset_index(drop=True)

    merged_df = pd.concat([meta, feats_df], axis=1)
    merged_df = merged_df[merged_df[eval_key].notna()].copy()
    feature_cols = feats_df.columns

    merged_df = merged_df.dropna(subset=feature_cols)  # drop NAs since we cannot compute
    asw = silhouette_score(merged_df[feature_cols].values, merged_df[eval_key].values, metric="cosine")
    asw = (asw + 1) / 2  # normalizes to [0, 1]

    np.save(asw_path, asw)



def silhouette_batch(parquet_path, label_key, batch_key, asw_batch_path):
    adata = filter_dmso_anndata(parquet_path)
    adata.obsm["X_embd"] = adata.X
    asw_batch = metrics.silhouette_batch(
        adata,
        batch_key,
        label_key,
        "X_embd",
        metric="cosine",
        verbose=False,
    )
    np.save(asw_batch_path, asw_batch)


def pcr_batch(pre_parquet_path, post_parquet_path, batch_key, pcr_batch_path):
    adata = filter_dmso_anndata(pre_parquet_path)
    adata_int = filter_dmso_anndata(post_parquet_path)
    pcr_score = metrics.pcr_comparison(
        adata,
        adata_int,
        embed=None,
        covariate=batch_key,
        verbose=False
    )

    np.save(pcr_batch_path, pcr_score)            


def pcr(parquet_path, batch_key, pcr_path):
    meta, vals, _ = filter_dmso(parquet_path)
    # 1 - pcr to make it 0 worst 1 best
    pcr_score = 1 - pc_regression(vals, meta[batch_key].values)
    np.save(pcr_path, pcr_score)


def isolated_labels_asw(adata_path, label_key, batch_key, il_asw_path):
    adata = ad.read_h5ad(adata_path)
    adata.obsm["X_embd"] = adata.X
    il_score_asw = metrics.isolated_labels(
        adata,
        label_key=label_key,
        batch_key=batch_key,
        embed="X_embd",
        cluster=False,
        # max number of batches a label should be present to be considered isolated
        iso_threshold=None,
        verbose=False,
    )
    np.save(il_asw_path, il_score_asw)


def isolated_labels_f1(adata_path, label_key, batch_key, il_f1_path):
    adata = ad.read_h5ad(adata_path)
    adata.obsm["X_embd"] = adata.X
    il_score_f1 = metrics.isolated_labels(
        adata,
        label_key=label_key,
        batch_key=batch_key,
        embed="X_embd",
        cluster=True,
        # max number of batches a label should be present to be considered isolated
        iso_threshold=None,
        verbose=False,
    )
    np.save(il_f1_path, il_score_f1)


def graph_connectivity(adata_path, label_key, graph_conn_path):
    adata = ad.read_h5ad(adata_path)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",
                                category=FutureWarning,
                                message="pandas.value_counts is deprecated")
        graph_conn_score = metrics.graph_connectivity(adata, label_key=label_key)
        np.save(graph_conn_path, graph_conn_score)


def kbet(adata_path, label_key, batch_key, kbet_path):
    adata = ad.read_h5ad(adata_path)
    M = metrics.kbet.diffusion_conn(adata, min_k=15, copy=False)
    adata.obsp["connectivities"] = M
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", module="rpy2.robjects")
        warnings.filterwarnings("ignore", category=FutureWarning)
        kbet_score = metrics.kBET(
            adata,
            batch_key=batch_key,
            label_key=label_key,
            type_="knn",
            embed=None,
            scaled=True,
            verbose=False,
        )
        np.save(kbet_path, kbet_score)


def lisi_label(adata_path, label_key, lisi_label_path):
    adata = ad.read_h5ad(adata_path)
    clisi = metrics.clisi_knn(
        adata,
        label_key=label_key,
        type_="knn",
        subsample=100,  # Use all data
        scale=True,
        n_cores=8,
        verbose=False,
    )
    np.save(lisi_label_path, clisi)


def lisi_batch(adata_path, batch_key, lisi_batch_path):
    adata = ad.read_h5ad(adata_path)
    
    ilisi = metrics.ilisi_knn(
        adata,
        batches=adata.obs[batch_key],
        scale=True,
    )
    print(ilisi)
    np.save(lisi_batch_path, ilisi)


def concat_metrics_per_method_into_df(metric_paths, scib_metrics, eval_keys, output_path):
    '''Concatenate scib metrics in a single file'''

    metrics_common_path = os.path.commonprefix(metric_paths)

    metrics = []
    for metric_path in metric_paths:
        # if the path contains an eval_key (f.e. 'asw_Metadata_MOA') we first extract the eval_key
        # then we remove it to end up with the metric. If we don't have an eval_key, the substring
        # is directly the metric name
        metric_substring = metric_path.replace(metrics_common_path, "").replace(".npy", "")
        eval_col = next((ek for ek in eval_keys if ek in metric_substring), None)
        if eval_col:
            metric = metric_substring.replace(f"_{eval_col}", "")
        else:
            metric = metric_substring
        score = np.load(metric_path) or None
        
        metrics.append((metric, eval_col, score))

    df = pd.DataFrame(metrics, columns=["metric", "eval_key", "score"])
    # df["eval_key"] = df["eval_key"].where(pd.notnull(df["eval_key"]), pd.NA)
    # df["score"] = df["score"].where(pd.notnull(df["score"]), pd.NA)
    df["score"] = df["score"].astype("float64")

    df.to_parquet(output_path)
