import pandas as pd
import numpy as np

from preprocessing import io


def mnn(parquet_path, batch_key, output_path):
    from mnnpy import mnn_correct, settings
    # disabling parallel because it silently freezes the execution
    settings.normalization = 'no_parallel'

    meta, vals, features = io.split_parquet(parquet_path)

    # group values
    col = pd.Series(data=meta.index, index=meta[batch_key])
    indices = list(col.groupby(level=0).indices.values())
    vals = [vals[ix] for ix in indices]

    # correct
    vals, _, __ = mnn_correct(*vals, var_index=features, do_concatenate=True)

    # Recover original order for vals
    indices = np.concatenate(indices)
    sort_ix = np.argsort(indices)
    vals = vals[sort_ix]

    # Save file
    features = [f'mnn_{i}' for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)
