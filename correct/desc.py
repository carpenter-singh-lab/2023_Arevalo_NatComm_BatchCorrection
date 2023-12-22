import scib
from quality_control import io


def desc(dframe_path: str, batch_key: str, output_path: str):
    '''DESC correction'''
    adata = io.to_anndata(dframe_path)
    meta = adata.obs.reset_index(drop=True).copy()

    adata_corrected = scib.integration.desc(adata, batch_key, use_gpu=True, gpu_id=0)
    vals = adata_corrected.obsm['X_emb']

    features = [f'desc_{i}' for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)
