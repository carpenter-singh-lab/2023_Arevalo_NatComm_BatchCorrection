import argparse
import logging
from scvi.external import SysVI
from preprocessing import io

logger = logging.getLogger(__name__)


def correct_with_sysvi(
    dframe_path: str,
    batch_key: list[str] | str,
    label_key: str,
    output_path: str,
    smoketest: bool = False,
):
    """sysVI correction from https://github.com/Hrovatin/scvi-tutorials/blob/main/scrna/sysVI.ipynb"""
    n_latent = 30
    n_epochs = 2 if smoketest else 999999

    adata = io.to_anndata(dframe_path)
    meta = adata.obs.reset_index(drop=True).copy()

    if isinstance(batch_key, list) and len(batch_key) == 1:
        batch_key = batch_key[0]

    batch_key = batch_key.split(",")

    if isinstance(batch_key, list):
        actual_batch_key = batch_key[0]
        assert isinstance(actual_batch_key, str)

        categorical_covariate_keys = batch_key[1:]
        if isinstance(categorical_covariate_keys, str):
            categorical_covariate_keys = [categorical_covariate_keys]
    else:
        actual_batch_key = batch_key
        categorical_covariate_keys = [None]

    SysVI.setup_anndata(
        adata,
        batch_key=actual_batch_key,
        categorical_covariate_keys=categorical_covariate_keys,
    )
    vae = SysVI(adata, n_layers=2, n_latent=n_latent, prior="standard_normal")
    vae.view_anndata_setup(adata=adata)
    vae.train(
        max_epochs=n_epochs,
        early_stopping=True,
        early_stopping_monitor="validation_loss",
    )

    vals = vae.get_latent_representation(adata=adata)
    features = [f"sysvi_{i}" for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform sysVI correction on data.")

    parser.add_argument("--input_data", required=True, help="Path to input data")
    parser.add_argument(
        "--batch_key",
        required=True,
        nargs="+",
        help="Batch key(s); provide multiple keys separated by spaces or as a comma-separated string.",
    )
    parser.add_argument("--label_key", required=True, help="Label key")
    parser.add_argument(
        "--output_path", required=True, help="Path to save corrected data"
    )
    parser.add_argument(
        "--smoketest",
        action="store_true",
        help="Run a smoketest with limited epochs (sets max_epochs to 2)",
    )

    args = parser.parse_args()

    correct_with_sysvi(
        dframe_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        output_path=args.output_path,
        smoketest=args.smoketest,
    )
