import logging
import argparse
from typing import List
from pycytominer.operations import Spherize
from preprocessing.io import merge_parquet, split_parquet

import numpy as np

logger = logging.getLogger(__name__)


def correct_with_sphering(
    dframe_path: str,
    method: str,
    column_norm: str,
    values_norm: List[str] | str,
    output_path: str,
    epsilon: float = 1e-6,
    **kwargs,
):
    spherer = Spherize(epsilon=epsilon, method=method)
    meta, vals, features = split_parquet(dframe_path)
    train_ix = meta[column_norm].isin(values_norm).values
    spherer.fit(vals[train_ix])
    vals = spherer.transform(vals).astype(np.float32)
    merge_parquet(meta, vals, features, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform sphering correction on data.")

    parser.add_argument("--input_data", required=True, help="Path to input data")
    parser.add_argument("--method", default="ZCA-cor", help="Sphering method")
    parser.add_argument(
        "--epsilon",
        type=float,
        default=0.1,
        help="Regularization parameter for sphering",
    )
    parser.add_argument(
        "--column_norm", required=True, help="Column name to use for normalization"
    )
    parser.add_argument(
        "--values_norm",
        required=True,
        help="Values in column_norm to use for normalization, separated by spaces",
    )
    parser.add_argument(
        "--output_path", required=True, help="Path to save sphered data"
    )

    args = parser.parse_args()

    if isinstance(args.values_norm, str):
        values_norm = [args.values_norm]

    correct_with_sphering(
        dframe_path=args.input_data,
        method=args.method,
        epsilon=args.epsilon,
        column_norm=args.column_norm,
        values_norm=values_norm,
        output_path=args.output_path,
    )
