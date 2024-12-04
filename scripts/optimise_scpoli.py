import os
import sys
import logging
import argparse
from preprocessing import io
import pandas as pd
import anndata as ad
import optuna
from scib_metrics.benchmark import Benchmarker
from scarches.models.scpoli import scPoli

logger = logging.getLogger(__name__)

def scib_benchmark_embedding(
    adata: ad.AnnData,
    batch_key: str,
    label_key: str,
) -> float:
    adata.obsm["trial"] = adata.X

    # silence output
    sys.stdout = open(os.devnull, "w")

    bm = Benchmarker(
        adata=adata,
        batch_key=batch_key,
        label_key=label_key,
        embedding_obsm_keys=["trial"],
    )
    bm.benchmark()
    df = bm.get_results(min_max_scale=False)

    # restore output
    sys.stdout.close()
    sys.stdout = sys.__stdout__

    return df.loc["trial"][["Batch correction", "Bio conservation"]].values


def objective(
    trial,
    adata: ad.AnnData, 
    batch_key: str, 
    label_key: str, 
    smoketest: bool = False,
):

    # silence output
    sys.stdout = open(os.devnull, "w")

    # Optimize hidden layer sizes
    num_layers = trial.suggest_int("num_layers", 1, 3)  # 1 to 4 layers
    hidden_layer_sizes = [trial.suggest_int(f"layer_{i}_size", 32, 512, step=32) for i in range(num_layers)]

    # Optimize latent dimensions and embedding size
    latent_dim = trial.suggest_int("latent_dim", 16, 128, step=16)
    embedding_dims = trial.suggest_int("embedding_dims", 2, 20, step=1)

    # Optimize pretraining to training epoch ratio
    total_epochs = 4 if smoketest else 50
    pretrain_to_train_ratio = trial.suggest_float("pretrain_to_train_ratio", 0.1, 0.9, step=0.1)
    n_pretrain_epochs = int(total_epochs * pretrain_to_train_ratio)
    n_train_epochs = total_epochs - n_pretrain_epochs

    # Optimize other parameters
    alpha_epoch_anneal = trial.suggest_int("alpha_epoch_anneal", 100, 1000, step=100)
    eta = trial.suggest_float("eta", 0.1, 1.0, step=0.1)

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

    integrated_adata = ad.AnnData(
        X=pd.DataFrame(vals, columns=features, index=adata.obs_names),
        obs=adata.obs,
    )

    batch, bio = scib_benchmark_embedding(
        adata=integrated_adata,
        batch_key=batch_key, 
        label_key=label_key
    )

    # restore output
    sys.stdout.close()
    sys.stdout = sys.__stdout__

    return batch, bio

def optimize_scpoli(
    input_path: str,
    batch_key: str,
    label_key: str,
    n_trials: int,
    output_path: str,
    smoketest: bool = False,
):
    if smoketest:
        n_trials = 2
    adata = io.to_anndata(input_path)

    study = optuna.create_study(directions=["maximize", "maximize"])
    study.optimize(lambda trial: objective(trial, adata.copy(), batch_key, label_key, smoketest), n_trials=n_trials)

    df = study.trials_dataframe()
    df = df.rename(columns={"values_0": "batch", "values_1": "bio"})
    df["total"] = 0.6 * df["bio"] + 0.4 * df["batch"]
    df = df.sort_values("total", ascending=False)
    df.to_csv(output_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Use Optuna to tune hyperparameters for scPoli.")

    parser.add_argument("--input_data", required=True, help="Path to input data.")
    parser.add_argument("--batch_key", required=True, help="Batch key.")
    parser.add_argument("--label_key", required=True, help="Label key.")
    parser.add_argument("--n_trials", required=True, help="How many trials to run.")
    parser.add_argument("--output_path", required=True, help="Where to save the optimal parameter set.")
    parser.add_argument(
        "--smoketest", action="store_true", help="Run a smoketest with limited epochs"
    )

    args = parser.parse_args()

    optimize_scpoli(
        input_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        n_trials=int(args.n_trials),
        output_path=args.output_path,
        smoketest=args.smoketest,
    )
