import logging
import argparse
from typing import List, Union
from scarches.models.scpoli import scPoli
from preprocessing import io

logger = logging.getLogger(__name__)


def correct_with_scpoli(
    dframe_path: str,
    batch_key: Union[List[str], str],  # we're using Python 3.9
    label_key: str,
    output_path: str,
    smoketest: bool = False,
    **kwargs,
):
    """scPoli correction from https://www.nature.com/articles/s41592-023-02035-2"""
    n_latent = 30
    n_epochs = (4, 2) if smoketest else (999999, 25)

    if isinstance(batch_key, list) and len(batch_key) == 1:
        if "," in batch_key[0]:
            batch_key = batch_key[0].split(",")

    if isinstance(batch_key, str) and "," in batch_key:
        batch_key = batch_key.split(",")
    elif isinstance(batch_key, str):
        batch_key = [batch_key]
    # batch_key = batch_key.split(",") # returns a list with len = 1 if no comma present, saves conversion

    adata = io.to_anndata(dframe_path)
    meta = adata.obs.reset_index(drop=True).copy()

    model = scPoli(
        adata=adata,
        condition_keys=batch_key,
        cell_type_keys=label_key,
        hidden_layer_sizes=[128],
        latent_dim=n_latent,
        embedding_dims=5,
        recon_loss="mse",
    )

    model.train(
        n_epochs=n_epochs[0],  # train until early stopping limit
        pretraining_epochs=n_epochs[1],
        use_early_stopping=True,
        alpha_epoch_anneal=1000,
        eta=0.5,
    )

    model.model.eval()
    vals = model.get_latent(adata, mean=True)
    features = [f"scpoli_{i}" for i in range(vals.shape[1])]
    print(meta, vals.shape, features, output_path)
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
    parser.add_argument(
        "--output_path", required=True, help="Path to save the output data"
    )
    parser.add_argument(
        "--smoketest", action="store_true", help="Run a smoketest with limited epochs"
    )

    args = parser.parse_args()

    correct_with_scpoli(
        dframe_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        output_path=args.output_path,
        smoketest=args.smoketest,
    )
