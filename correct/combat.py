import scanpy as sc
from preprocessing import io


def combat(dframe_path: str, batch_key: str, output_path: str):
    '''Combat correction'''

    adata = io.to_anndata(dframe_path)
    vals = sc.pp.combat(adata, key=batch_key, inplace=False)

    meta = adata.obs.reset_index(drop=True).copy()
    features = [f'combat_{i}' for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)
