import sys
import logging
from typing import Union, List
from scarches.models.scpoli import scPoli 
from preprocessing import io

logger = logging.getLogger(__name__)

def correct_with_scpoli(dframe_path: str, batch_key: Union[List[str], str], label_key: str, output_path: str, **kwargs):
    '''scPoli correction from https://www.nature.com/articles/s41592-023-02035-2'''
    n_latent = 30
    smoketest = kwargs.get("smoketest", 0)
    n_epochs = (999999, 25) if smoketest else (4, 2)

    batch_key = batch_key.split(' ') # returns a list with len = 1 if no whitespace present, saves conversion

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
        n_epochs=n_epochs[0], # train until early stopping limit
        pretraining_epochs=n_epochs[1],
        use_early_stopping=True,
        alpha_epoch_anneal=1000,
        eta=0.5,
    ) 

    model.model.eval()
    vals = model.get_latent(adata, mean=True)
    features = [f'scpoli_{i}' for i in range(vals.shape[1])]
    print(meta, vals.shape, features, output_path)
    io.merge_parquet(meta, vals, features, output_path)

if __name__ == "__main__":
    dframe_path = sys.argv[1]
    batch_key = sys.argv[2]
    label_key = sys.argv[3]
    output_path = sys.argv[4]
    correct_with_scpoli(dframe_path, batch_key, label_key, output_path)
