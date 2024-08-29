import sys
import logging
from desc import train, scale_bygroup
from preprocessing import io
import tempfile

logger = logging.getLogger(__name__)

def correct_with_desc(dframe_path: str, batch_key: str, output_path: str):
    '''DESC correction'''

    adata = io.to_anndata(dframe_path)
    meta = adata.obs.reset_index(drop=True).copy()

    res = 1.0
    scale_bygroup(adata, batch_key, max_value=None)
    adata = train(
        adata,
        dims=[adata.shape[1], 128, 32],
        tol=0.05,
        n_neighbors=64,
        batch_size=1024,
        louvain_resolution=res,
        save_encoder_weights=False,
        save_dir=tempfile.TemporaryDirectory().name,
        do_tsne=False,
        use_GPU=True,
        GPU_id=0,
        num_Cores=100,
        use_ae_weights=False,
        do_umap=False,
    )

    vals = adata.obsm[f"X_Embeded_z{res}"]

    features = [f'desc_{i}' for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)


if __name__ == "__main__":
    dframe_path = sys.argv[1]
    batch_key = sys.argv[2]
    output_path = sys.argv[3]

    correct_with_desc(dframe_path, batch_key, output_path)
