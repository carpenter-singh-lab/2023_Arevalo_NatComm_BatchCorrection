import sys
import logging
import pandas as pd
import numpy as np
import scanpy as sc
from preprocessing import io

logger = logging.getLogger(__name__)

def correct_with_scanorama(parquet_path, batch_key, output_path):
    '''Scanorama correction'''
    adata = io.to_anndata(parquet_path)
    # Scanorama requires samples sorted per batch_key.
    idx = adata.obs[batch_key].argsort().values
    adata = adata[idx, :]
    n_latent = min(adata.shape) - 1  # required for arpack
    sc.tl.pca(adata, n_comps=n_latent)  # required for arpack
    sc.external.pp.scanorama_integrate(adata,
                                       batch_key,
                                       adjusted_basis='X_scanorama')
    meta = adata.obs.reset_index(drop=True).copy()
    vals = adata.obsm['X_scanorama']
    features = [f'scanorama_{i}' for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)

def correct_with_pca_scanorama(parquet_path, batch_key, output_path):
    '''Scanorama correction'''
    adata = io.to_anndata(parquet_path)
    # Scanorama requires samples sorted per batch_key.
    idx = adata.obs[batch_key].argsort().values
    adata = adata[idx, :]
    n_latent = min(adata.shape) - 1  # required for arpack
    sc.tl.pca(adata, n_comps=n_latent)  # required for arpack
    sc.external.pp.scanorama_integrate(adata,
                                       batch_key,
                                       adjusted_basis='X_scanorama')
    meta = adata.obs.reset_index(drop=True).copy()
    vals = adata.obsm['X_scanorama']
    features = [f'scanorama_{i}' for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)



if __name__ == "__main__":
    mode = sys.argv[1]
    dframe_path = sys.argv[2]
    batch_key = sys.argv[3]
    output_path = sys.argv[4]

    if mode == "scanorama":
        correct_with_scanorama(dframe_path, batch_key, output_path)
    elif mode == "pca_scanorama":
        correct_with_pca_scanorama(dframe_path, batch_key, output_path)
    else:
        raise ValueError("Invalid mode. Choose either 'scanorama' or 'pca_scanorama'")
