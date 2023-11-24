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
sc.pp.neighbors(cl_adata, use_rep='X', n_neighbors=15)
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

pcr_score = metrics.pcr_comparison(
    raw_adata, adata, embed=None, covariate=batch_key, verbose=False
)
