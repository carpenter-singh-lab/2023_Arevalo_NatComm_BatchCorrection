import logging
import argparse
import scanpy as sc
from preprocessing import io

logger = logging.getLogger(__name__)


def correct_with_combat(dframe_path: str, batch_key: str, output_path: str):
    """Combat correction"""

    adata = io.to_anndata(dframe_path)
    vals = sc.pp.combat(adata, key=batch_key, inplace=False)

    meta = adata.obs.reset_index(drop=True).copy()
    features = [f"combat_{i}" for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)        
    logger.info(f"Combat-corrected data saved successfully to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform Combat correction on data.")

    parser.add_argument("--input_data", required=True, help="Path to input data.")
    parser.add_argument("--batch_key", required=True, help="Batch key.")
    parser.add_argument(
        "--output_path", required=True, help="Path to save corrected data."
    )

    args = parser.parse_args()

    correct_with_combat(
        dframe_path=args.input_data,
        batch_key=args.batch_key,
        output_path=args.output_path,
    )
