from scvi.external import SysVI
from preprocessing import io

def sysvi(dframe_path: str, batch_key: list[str] | str, label_key: str, output_path: str):
    '''sysVI correction from https://github.com/Hrovatin/scvi-tutorials/blob/main/scrna/sysVI.ipynb'''
    n_latent = 30

    adata = io.to_anndata(dframe_path)
    meta = adata.obs.reset_index(drop=True).copy()

    # min_value = adata.X.min()
    # adata.X -= min_value

    if isinstance(batch_key, list):
        actual_batch_key = batch_key[0]
    else:
        categorical_covariate_keys = batch_key[1:]

    SysVI.setup_anndata(adata, batch_key=batch_key, categorical_covariate_keys=categorical_covariate_keys)
    vae = SysVI(adata, n_layers=2, n_latent=n_latent, prior="standard_normal")
    vae.view_anndata_setup(adata=adata)
    vae.train()

    vals = vae.get_latent_representation()
    features = [f'sysvi_{i}' for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)
