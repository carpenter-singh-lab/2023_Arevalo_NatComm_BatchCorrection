import anndata as ad
import scanpy as sc
from bbknn.matrix import bbknn

from preprocessing import io


def clustering(input_path, output_path, batch_key):
    """Compute metrics for BBKNN algorithm. Separate function because not all
    metrics used in the benchmark can be computed"""
    meta, vals, features = io.split_parquet(input_path)
    meta.index = meta.index.astype(str)
    neighbors_within_batch = 3 if len(vals) < 1e5 else 25
    bbknn_out = bbknn(
        vals, meta[batch_key].values, neighbors_within_batch=neighbors_within_batch
    )

    adata = ad.AnnData(vals, obs=meta)
    adata.obsp["distances"] = bbknn_out[0]
    adata.obsp["connectivities"] = bbknn_out[1]
    adata.uns["neighbors"] = {
        "distances_key": "distances",
        "connectivities_key": "connectivities",
    }
    adata.X = None
    sc.tl.leiden(adata, key_added="Metadata_Cluster")
    adata.write_h5ad(output_path, compression="gzip")
