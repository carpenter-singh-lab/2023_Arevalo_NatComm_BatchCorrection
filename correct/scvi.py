from scvi.model import SCVI
from preprocessing import io


def scvi(dframe_path: str, batch_key: str, label_key: str, output_path: str):
    '''scVI correction'''
    n_latent = 30

    adata = io.to_anndata(dframe_path)
    meta = adata.obs.reset_index(drop=True).copy()

    min_value = adata.X.min()
    adata.X -= min_value

    SCVI.setup_anndata(adata, batch_key=batch_key, labels_key=label_key)
    vae = SCVI(adata, n_layers=2, n_latent=n_latent)
    vae.view_anndata_setup(adata=adata)
    vae.train()

    vals = vae.get_latent_representation()
    features = [f'scvi_{i}' for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)
