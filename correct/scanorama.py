import pandas as pd
import numpy as np
import scanpy as sc
from preprocessing import io


def pca_scanorama(parquet_path, batch_key, output_path):
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


def scanorama(parquet_path, batch_key, output_path):
    '''Scanorama correction without pca'''
    from scanorama import assemble
    meta, vals, features = io.split_parquet(parquet_path)

    # group values
    col = pd.Series(data=meta.index, index=meta[batch_key])
    indices = list(col.groupby(level=0).indices.values())
    vals = [vals[ix] for ix in indices]

    # correct
    vals = assemble(vals, batch_size=5000)
    vals = np.concatenate(vals)

    # Recover original order for vals
    indices = np.concatenate(indices)
    sort_ix = np.argsort(indices)
    vals = vals[sort_ix]

    # Save file
    features = [f'scanorama_{i}' for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)
