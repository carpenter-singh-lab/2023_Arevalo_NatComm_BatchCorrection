import warnings

import anndata as ad
import scanpy as sc

from quality_control import io


def mde(parquet_path, mde_path):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=UserWarning,
            message="Failed to load image Python extension",
        )
        import pymde
    meta, vals, _ = io.split_parquet(parquet_path)
    # params taken from scvi.model.utils.mde
    mde_params = {
        "embedding_dim": 2,
        "constraint": pymde.Standardized(),
        "repulsive_fraction": 0.7,
        "verbose": False,
        "n_neighbors": 15,
    }

    embd = pymde.preserve_neighbors(vals, **mde_params).embed()

    meta["x"] = embd[:, 0]
    meta["y"] = embd[:, 1]
    meta.to_parquet(mde_path)


def pca(parquet_path, pca_path):
    meta, vals, _ = io.split_parquet(parquet_path)
    embd = sc.tl.pca(vals, n_comps=2)  # Generates X_pca
    meta["x"] = embd[:, 0]
    meta["y"] = embd[:, 1]
    meta.to_parquet(pca_path)


def umap(adata_path, umap_path):
    adata = ad.read_h5ad(adata_path)
    sc.tl.umap(adata)  # Generates X_umap
    meta = adata.obs
    meta.reset_index(drop=True, inplace=True)
    meta["x"] = adata.obsm["X_umap"][:, 0]
    meta["y"] = adata.obsm["X_umap"][:, 1]
    meta.to_parquet(umap_path)
