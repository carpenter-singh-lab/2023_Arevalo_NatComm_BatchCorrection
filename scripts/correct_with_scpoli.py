import logging
import argparse
from typing import List, Union, Optional, Literal
from scarches.models.scpoli import scPoli
from preprocessing import io
import scanpy as sc
import anndata as ad
import pandas as pd

logger = logging.getLogger(__name__)


def correct_with_scpoli(
    dframe_path: str,
    batch_key: Union[List[str], str],  # we're using Python 3.9
    label_key: str,
    parameter_path: str,
    output_path: str,
    preproc: Optional[Literal["pca"]] = None,
    smoketest: bool = False,
    **kwargs,
):
    """scPoli correction from https://www.nature.com/articles/s41592-023-02035-2"""

    # load hyperparameters
    params = pd.read_csv(parameter_path)
    params = params.sort_values("total").iloc[0].to_dict()

    if params["state"] != "COMPLETE":
        raise ValueError("Optimization did not complete successfully")

    alpha_epoch_anneal = params["params_alpha_epoch_anneal"]
    embedding_dims = params["params_embedding_dims"]
    eta = params["params_eta"]
    latent_dim = params["params_latent_dim"]
    hidden_layer_sizes = [
        params[f"params_layer_{i}_size"] for i in range(params["params_num_layers"])
    ]
    pretrain_to_train_ratio = params["params_pretrain_to_train_ratio"]

    total_epochs = 4 if smoketest else 100
    n_pretrain_epochs = 2 if smoketest else int(total_epochs * pretrain_to_train_ratio)
    n_train_epochs = 2 if smoketest else (total_epochs - n_pretrain_epochs)
    n_train_epochs += (0 if smoketest else 9999999) # we train until convergence

    print("\nUsing the following hyperparameters:")
    print(f"- alpha_epoch_anneal: {alpha_epoch_anneal}")
    print(f"- embedding_dims: {embedding_dims}")
    print(f"- eta: {eta}")
    print(f"- latent_dim: {latent_dim}")
    print(f"- hidden_layer_sizes: {hidden_layer_sizes}")
    print(f"- n_pretrain_epochs: {n_pretrain_epochs}")
    print(f"- n_train_epochs: {n_train_epochs}\n")

    if isinstance(batch_key, list) and len(batch_key) == 1 and "," in batch_key[0]:
        batch_key = batch_key[0].split(",")

    if isinstance(batch_key, str):
        batch_key = batch_key.split(",") if "," in batch_key else [batch_key]

    adata = io.to_anndata(dframe_path)
    meta = adata.obs.reset_index(drop=True).copy()

    if preproc == "pca":
        logger.info("Applying PCA preprocessing with Scanpy")
        sc.pp.pca(adata, svd_solver="arpack", n_comps=50)
        # Create new adata object with PCA-transformed data
        adata_for_training = ad.AnnData(adata.obsm["X_pca"].copy())
        adata_for_training.obs = adata.obs.copy()
        adata = adata_for_training
        logger.info("PCA completed. Shape of PCA-transformed data: %s", adata.X.shape)


    model = scPoli(
        adata=adata,
        condition_keys=batch_key,
        cell_type_keys=label_key,
        hidden_layer_sizes=hidden_layer_sizes,
        latent_dim=latent_dim,
        embedding_dims=embedding_dims,
        recon_loss="mse",
    )

    model.train(
        n_epochs=n_train_epochs,
        pretraining_epochs=n_pretrain_epochs,
        use_early_stopping=True,
        alpha_epoch_anneal=alpha_epoch_anneal,
        eta=eta,
    )

    model.model.eval()
    vals = model.get_latent(adata, mean=True)
    features = [f"scpoli_{i}" for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="scPoli correction from https://www.nature.com/articles/s41592-023-02035-2"
    )

    parser.add_argument("--input_data", required=True, help="Path to the input data")
    parser.add_argument(
        "--batch_key",
        nargs="+",
        required=True,
        help="Batch key(s), provide multiple keys separated by commas",
    )
    parser.add_argument("--label_key", required=True, help="Label key")
    parser.add_argument("--parameter_path", required=True, help="Path to the parameter file")
    parser.add_argument(
        "--output_path", required=True, help="Path to save the output data"
    )
    parser.add_argument(
        "--preproc",
        type=str,
        choices=["pca"],
        default=None,
        help="Preprocessing method to apply. Choices: 'pca' or None",
    )
    parser.add_argument(
        "--smoketest", action="store_true", help="Run a smoketest with limited epochs"
    )

    args = parser.parse_args()

    correct_with_scpoli(
        dframe_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        parameter_path=args.parameter_path,
        output_path=args.output_path,
        preproc=args.preproc,
        smoketest=args.smoketest,
    )
