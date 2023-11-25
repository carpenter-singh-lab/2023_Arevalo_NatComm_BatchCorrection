import numpy as np
import anndata as ad
import scanpy as sc
from scvi.model.utils import mde as scvi_mde
from quality_control import io


def mde(parquet_path, mde_path):
    meta, vals, _ = io.split_parquet(parquet_path)
    embd = scvi_mde(vals)
    meta['x'] = embd[:, 0]
    meta['y'] = embd[:, 1]
    meta.to_parquet(mde_path)


def pca(parquet_path, pca_path):
    meta, vals, _ = io.split_parquet(parquet_path)
    embd = sc.tl.pca(vals, n_comps=2)  # Generates X_pca
    meta['x'] = embd[:, 0]
    meta['y'] = embd[:, 1]
    meta.to_parquet(pca_path)


def umap(parquet_path, umap_path):
    # TODO: to avoid neighbors computation use instead in the smk rule
    # input: f"outputs/{{scenario}}/metrics/{criteria}/scib/{{pipeline}}_clusters.h5ad",
    adata = io.to_anndata(parquet_path)
    sc.pp.neighbors(adata, use_rep='X')
    sc.tl.umap(adata)  # Generates X_umap
    import ipdb; ipdb.set_trace() # BREAKPOINT
    meta = adata.obs
    meta.reset_index(drop=True, inplace=True)
    meta['x'] = adata.obsm['X_umap'][:, 0]
    meta['y'] = adata.obsm['X_umap'][:, 1]
    meta.to_parquet(umap_path)
