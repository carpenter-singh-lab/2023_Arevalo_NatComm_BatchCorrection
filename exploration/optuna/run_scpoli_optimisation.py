import os
import sys
import pandas as pd
import anndata as ad
import optuna
import pickle
from scarches.models.scpoli import scPoli
from scib_metrics.benchmark import Benchmarker

data_path = "/home/icb/tim.treis/projects/broad_integrate/2023_Arevalo_BatchCorrection/outputs/scenario_7/mad_int_featselect.parquet"

data = pd.read_parquet(data_path)
metadata_cols = data.filter(regex="Metadata").columns

adata = ad.AnnData(X=data.drop(metadata_cols, axis=1).values, obs=data[metadata_cols])
adata

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
        n_jobs=4,
    )
    bm.benchmark()
    df = bm.get_results(min_max_scale=False)

    # restore output
    sys.stdout.close()
    sys.stdout = sys.__stdout__

    return df.loc["trial"][["Batch correction", "Bio conservation"]].values


def objective(trial, adata, batch_key, label_key):

    # silence output
    sys.stdout = open(os.devnull, "w")

    # Optimize hidden layer sizes
    num_layers = trial.suggest_int("num_layers", 1, 3)  # 1 to 4 layers
    hidden_layer_sizes = [trial.suggest_int(f"layer_{i}_size", 32, 512, step=32) for i in range(num_layers)]
    
    # Optimize latent dimensions and embedding size
    latent_dim = trial.suggest_int("latent_dim", 16, 128, step=16)
    embedding_dims = trial.suggest_int("embedding_dims", 2, 20, step=1)

    # Optimize pretraining to training epoch ratio
    total_epochs = 75
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
        label_key=label_key,
    )

    # restore output
    sys.stdout.close()
    sys.stdout = sys.__stdout__

    return batch, bio

batch_key = "Metadata_Source"
label_key = "Metadata_JCP2022"

study = optuna.create_study(directions=["maximize", "maximize"])
study.optimize(lambda trial: objective(trial, adata.copy(), batch_key, label_key), n_trials=50)

with open("./optuna_study.pkl", "wb") as output_file:
    pickle.dump(study, output_file)