import numpy as np
from quality_control import io
from copairs.matching import reverse_index


def mnn(parquet_path, batch_key, output_path):
    from mnnpy import mnn_correct
    meta, vals, features = io.split_parquet(parquet_path)
    indices = reverse_index(meta[batch_key])
    vals = [vals[ix] for ix in indices.values]
    vals, _, __ = mnn_correct(*vals, var_index=features, do_concatenate=True)

    indices = np.concatenate(indices.values)
    meta = meta.iloc[indices]
    features = [f'mnn_{i}' for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)
