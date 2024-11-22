import argparse
import logging
import pandas as pd
import numpy as np
from typing import List, Union
from mnnpy import mnn_correct, settings
from preprocessing import io

logger = logging.getLogger(__name__)


def correct_with_mnn(
    parquet_path: str,
    batch_key: Union[str, List[str]],
    output_path: str,
):
    """MNN correction"""
    # disabling parallel because it silently freezes the execution
    settings.normalization = "no_parallel"

    meta, vals, features = io.split_parquet(parquet_path)
    vals = vals.astype(np.float32)

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
    features = [f"mnn_{i}" for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform MNN correction on data.")

    parser.add_argument(
        "--input_data", required=True, help="Path to the input data in Parquet format."
    )
    parser.add_argument("--batch_key", required=True, help="Batch key column name.")
    parser.add_argument(
        "--output_path", required=True, help="Path to save the corrected data."
    )

    args = parser.parse_args()

    correct_with_mnn(
        parquet_path=args.input_data,
        batch_key=args.batch_key,
        output_path=args.output_path,
    )
