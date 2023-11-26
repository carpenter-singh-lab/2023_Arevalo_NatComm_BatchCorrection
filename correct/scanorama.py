from scanorama import assemble
import numpy as np
import scanpy as sc
from quality_control import io
from copairs.matching import reverse_index


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
    meta, vals, features = io.split_parquet(parquet_path)
    indices = reverse_index(meta[batch_key])
    vals = [vals[ix] for ix in indices.values]
    vals = assemble(vals)
    vals = np.concatenate(vals)
    indices = np.concatenate(indices.values)
    meta = meta.iloc[indices]
    io.merge_parquet(meta, vals, features, output_path)
