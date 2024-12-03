import logging
import argparse
import scvi
from preprocessing import io

logger = logging.getLogger(__name__)

def correct_with_gaushanvi(
    dframe_path: str,
    batch_key: str | list,
    label_key: str,
    output_path: str,
    smoketest: bool = False,
):
    """gaushanvi correction"""
    n_latent = 30
    n_epochs = 2 if smoketest else 999999

    adata = io.to_anndata(dframe_path)
    meta = adata.obs.reset_index(drop=True).copy()

    # Adjust data (ensuring all values are non-negative)
    min_value = adata.X.min()
    adata.X -= min_value

    scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key, labels_key=label_key)
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=n_latent, gene_likelihood="normal")
    vae.train(
        max_epochs=n_epochs,
        early_stopping=True,
        early_stopping_monitor="elbo_validation",
    )

    # gashanvi setup
    scanvi_model = scvi.model.SCANVI.from_scvi_model(vae, unlabeled_category="unknown")
    scanvi_model.train(max_epochs=n_epochs)

    # Get latent representation
    vals = scanvi_model.get_latent_representation()
    features = [f'gaushanvi_{i}' for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform gaushanvi correction on data.")

    parser.add_argument("--input_data", required=True, help="Path to input data")
    parser.add_argument("--batch_key", required=True, help="Batch key")
    parser.add_argument("--label_key", required=True, help="Label key")
    parser.add_argument("--output_path", required=True, help="Path to save corrected data")
    parser.add_argument(
        "--smoketest",
        action="store_true",
        help="Run a smoketest with limited epochs",
    )

    args = parser.parse_args()

    correct_with_gaushanvi(
        dframe_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        output_path=args.output_path,
        smoketest=args.smoketest,
    )
