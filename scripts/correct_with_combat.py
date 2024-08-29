import sys
import logging
import scanpy as sc
from preprocessing import io

logger = logging.getLogger(__name__)

def correct_with_combat(dframe_path: str, batch_key: str, output_path: str):
    '''Combat correction'''

    adata = io.to_anndata(dframe_path)
    vals = sc.pp.combat(adata, key=batch_key, inplace=False)

    meta = adata.obs.reset_index(drop=True).copy()
    features = [f'combat_{i}' for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)


if __name__ == "__main__":
    dframe_path = sys.argv[1]
    batch_key = sys.argv[2]
    output_path = sys.argv[3]

    correct_with_combat(dframe_path, batch_key, output_path)
