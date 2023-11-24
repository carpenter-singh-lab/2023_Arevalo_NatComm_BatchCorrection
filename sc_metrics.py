import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

from quality_control import io
from scib import metrics
import scanpy as sc
import anndata as ad
from sklearn.metrics import silhouette_score

parquet_path = 'outputs/scenario_1/mad_featselect_sphering.parquet'
raw_parquet_path = 'outputs/scenario_1/mad_featselect.parquet'

label_key = 'Metadata_JCP2022'
batch_key = 'Metadata_Batch'
cluster_key = 'Metadata_Cluster'
meta, feats, features = io.split_parquet(parquet_path)
meta.index = meta.index.astype(str)
adata = ad.AnnData(feats, meta)
raw_meta, raw_feats, raw_features = io.split_parquet(raw_parquet_path)
raw_meta.index = raw_meta.index.astype(str)
raw_adata = ad.AnnData(raw_feats, raw_meta)

# Get cluster and neighbors
cl_adata = ad.AnnData(feats, meta)
print('neighbors')
sc.pp.neighbors(cl_adata, use_rep='X', n_neighbors=15, metric='cosine')
print('clustering')
metrics.cluster_optimal_resolution(cl_adata,
                                   label_key=label_key,
                                   cluster_key=cluster_key,
                                   metric=metrics.nmi,
                                   verbose=False)
adata.obs[cluster_key] = cl_adata.obs[cluster_key]

nmi = metrics.nmi(adata, label_key, cluster_key)
ari = metrics.ari(adata, label_key, cluster_key)

asw = silhouette_score(feats, meta[label_key], metric='cosine')
asw = (asw + 1) / 2

adata.obsm['X_sphe'] = feats
asw_batch = metrics.silhouette_batch(
    adata,
    batch_key,
    label_key,
    'X_sphe',
    metric="cosine",
    verbose=False,
)

pcr_score = metrics.pcr_comparison(raw_adata,
                                   adata,
                                   embed=None,
                                   covariate=batch_key,
                                   verbose=False)

il_score_asw = metrics.isolated_labels(
    adata,
    label_key=label_key,
    batch_key=batch_key,
    embed='X_sphe',
    cluster=False,
    iso_threshold=None, # max number of batches a label should be present to be considered isolated
    verbose=False,
)

il_score_f1 = metrics.isolated_labels(
    adata,
    label_key=label_key,
    batch_key=batch_key,
    embed='X_sphe',
    cluster=True,
    iso_threshold=None, # max number of batches a label should be present to be considered isolated
    verbose=True,
)
graph_conn_score = metrics.graph_connectivity(cl_adata, label_key=label_key)
import warnings
warnings.filterwarnings('ignore', module='rpy2.robjects')
warnings.filterwarnings('ignore', category=FutureWarning)
kbet_score = metrics.kBET(
            cl_adata,
            batch_key=batch_key,
            label_key=label_key,
            type_='knn',
            embed=None,
            scaled=True,
            verbose=False,
)
subsample = 1
clisi = metrics.clisi_graph(
    cl_adata,
    label_key=label_key,
    type_='knn',
    subsample=subsample * 100,
    scale=True,
    n_cores=8,
    verbose=True,
)

ilisi = metrics.ilisi_graph(
    cl_adata,
    batch_key=batch_key,
    type_='knn',
    subsample=subsample * 100,
    scale=True,
    n_cores=8,
    verbose=True,
)
