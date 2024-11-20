import argparse
import logging
import scanpy as sc
from preprocessing import io
import pandas as pd
import numpy as np
import anndata as ad

logger = logging.getLogger(__name__)


def correct_with_scanorama(parquet_path, batch_key, output_path, smoketest=False):
    """Scanorama correction on raw data"""
    import scanorama
    adata = io.to_anndata(parquet_path)

    def split_adata_by_col(adata, col):
        """Split the AnnData object by a column"""
        splits = []
        for value in adata.obs[col].cat.categories:
            mask = adata.obs[col] == value
            splits.append(adata[mask].copy())
        return splits

    adatas_by_source = split_adata_by_col(adata, batch_key)
    scanorama.integrate_scanpy(adatas_by_source)
    corrected_adata = ad.concat(adatas_by_source)

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
    features = [f"scanorama_{i}" for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)

def correct_with_pca_scanorama(parquet_path, batch_key, output_path, smoketest=False):
    """Scanorama correction using PCA components"""
    adata = io.to_anndata(parquet_path)

    # Sort adata based on batch_key (needed for scanorama)
    adata = adata[adata.obs.sort_values(by=batch_key).index]

    # Determine number of principal components
    n_pcs = 10 if smoketest else 50

    # Compute PCA
    sc.tl.pca(adata, n_comps=n_pcs)

    # Use the PCA components in Scanorama
    sc.external.pp.scanorama_integrate(
        adata, key=batch_key, dimred="X_pca", adjusted_basis="X_scanorama"
    )

    meta = adata.obs.reset_index(drop=True).copy()
    vals = adata.obsm["X_scanorama"]
    features = [f"scanorama_{i}" for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perform Scanorama correction on data."
    )

    parser.add_argument(
        "--mode",
        required=True,
        choices=["scanorama", "pca_scanorama"],
        help="Correction mode to use: 'scanorama' or 'pca_scanorama'.",
    )
    parser.add_argument(
        "--input_data", required=True, help="Path to the input data in Parquet format."
    )
    parser.add_argument("--batch_key", required=True, help="Batch key.")
    parser.add_argument(
        "--output_path", required=True, help="Path to save the corrected data."
    )
    parser.add_argument(
        "--smoketest",
        action="store_true",
        help="Run a smoketest with reduced computation for testing purposes.",
    )

    args = parser.parse_args()

    if args.mode == "scanorama":
        correct_with_scanorama(
            parquet_path=args.input_data,
            batch_key=args.batch_key,
            output_path=args.output_path,
            smoketest=args.smoketest,
        )
    elif args.mode == "pca_scanorama":
        correct_with_pca_scanorama(
            parquet_path=args.input_data,
            batch_key=args.batch_key,
            output_path=args.output_path,
            smoketest=args.smoketest,
        )
    else:
        raise ValueError(
            "Invalid mode. Choose either 'scanorama' or 'pca_scanorama'."
        )
