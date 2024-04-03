import anndata as ad
import pandas as pd
import scanpy as sc
from bbknn.matrix import bbknn
from scib import metrics

from metrics import DIMENSION_MAP
from preprocessing import io


def bbknn_metrics(input_path, output_path, batch_key):
    """Compute metrics for BBKNN algorithm. Separate function because not all
    metrics used in the benchmark can be computed"""
    meta, vals, features = io.split_parquet(input_path)
    meta.index = meta.index.astype(str)
    neighbors_within_batch = 3 if len(vals) < 1e5 else 25
    bbknn_out = bbknn(vals,
                      meta[batch_key].values,
                      neighbors_within_batch=neighbors_within_batch)

    adata = ad.AnnData(vals, obs=meta)
    adata.obsp["distances"] = bbknn_out[0]
    adata.obsp["connectivities"] = bbknn_out[1]
    adata.uns["neighbors"] = {
        "distances_key": "distances",
        "connectivities_key": "connectivities",
    }
    adata.X = None

    sc.tl.leiden(adata, key_added="Metadata_Cluster")
    scores = {}
    scores["graph_conn"] = metrics.graph_connectivity(adata,
                                                      "Metadata_JCP2022")
    scores["kbet"] = metrics.kbet.kBET(adata,
                                       "Metadata_Batch",
                                       "Metadata_JCP2022",
                                       type_="knn")
    scores["lisi_label"] = metrics.clisi_graph(adata,
                                               "Metadata_JCP2022",
                                               type_="knn")
    scores["lisi_batch"] = metrics.ilisi_graph(adata,
                                               "Metadata_Batch",
                                               type_="knn")
    scores["ari"] = metrics.ari(adata,
                                cluster_key="Metadata_Cluster",
                                label_key="Metadata_JCP2022")
    scores["nmi"] = metrics.nmi(adata,
                                cluster_key="Metadata_Cluster",
                                label_key="Metadata_JCP2022")

    scores = pd.Series(scores, name="score")
    scores.index.name = "metric"
    scores = scores.reset_index()
    scores["dimension"] = scores["metric"].map(DIMENSION_MAP)
    scores.to_parquet(output_path)
