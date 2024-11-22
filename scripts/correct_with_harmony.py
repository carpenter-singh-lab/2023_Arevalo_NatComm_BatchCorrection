import argparse
import logging
from harmonypy import run_harmony
from preprocessing import io
import scanpy as sc

logger = logging.getLogger(__name__)


def correct_with_harmony(
    dframe_path: str, batch_key: list, output_path: str, smoketest=False
):
    """Harmony correction"""
    n_max_iter = 2 if smoketest else 999999

    meta, feats, features = io.split_parquet(dframe_path)
    harmony_out = run_harmony(
        feats,
        meta,
        batch_key,
        max_iter_harmony=n_max_iter,
        nclust=300,  # Number of compounds
    )

    feats = harmony_out.Z_corr.T
    features = [f"harmony_{i}" for i in range(feats.shape[1])]
    io.merge_parquet(meta, feats, features, output_path)


def correct_with_harmony_pca(
    dframe_path: str, batch_key: list, output_path: str, smoketest=False
):
    """Harmony correction with PCA"""
    n_max_iter = 2 if smoketest else 999999

    meta, feats, features = io.split_parquet(dframe_path)
    n_latent = min(feats.shape) - 1  # Required for arpack
    logger.info("Computing PCA...")
    sc.tl.pca(feats, n_comps=n_latent)  # Generates X_pca
    logger.info("Computing PCA Done.")
    harmony_out = run_harmony(
        feats,
        meta,
        batch_key,
        max_iter_harmony=n_max_iter,
        nclust=300,  # Number of compounds
    )

    feats = harmony_out.Z_corr.T
    features = [f"harmony_{i}" for i in range(feats.shape[1])]
    io.merge_parquet(meta, feats, features, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform Harmony correction on data.")
    parser.add_argument(
        "--mode",
        choices=["harmony", "harmony_pca"],
        required=True,
        help="Correction mode to use.",
    )
    parser.add_argument("--input_data", required=True, help="Path to input data.")
    parser.add_argument("--batch_key", required=True, help="Batch key.")
    parser.add_argument(
        "--output_path", required=True, help="Path to save corrected data."
    )
    parser.add_argument(
        "--smoketest",
        action="store_true",
        help="Run in smoketest mode (limited iterations).",
    )

    args = parser.parse_args()

    if args.mode == "harmony":
        correct_with_harmony(
            dframe_path=args.input_data,
            batch_key=args.batch_key,
            output_path=args.output_path,
            smoketest=args.smoketest,
        )
    elif args.mode == "harmony_pca":
        correct_with_harmony_pca(
            dframe_path=args.input_data,
            batch_key=args.batch_key,
            output_path=args.output_path,
            smoketest=args.smoketest,
        )
    else:
        raise ValueError("Invalid mode. Choose either 'harmony' or 'harmony_pca'")
