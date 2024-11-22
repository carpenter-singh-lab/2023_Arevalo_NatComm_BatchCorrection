import sys
import logging
import scanpy as sc
import pandas as pd

from scib_metrics.benchmark import Benchmarker
from metrics.scib import _load_opentargets_moa_info, _merge_with_duplication

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

logger = logging.getLogger(__name__)


def run_scibmetrics_benchmarker(
    adata_path, output_paths, batch_key, eval_keys, methods
):
    adata = sc.read_h5ad(adata_path)

    if isinstance(batch_key, list) and len(batch_key) == 1:
        batch_key = batch_key[0]

    # keys are whitespace separated
    eval_keys = eval_keys.split(" ")

    results = []

    # iterate over keys and benchmark embeddings
    for eval_key in eval_keys:
        if eval_key == "Metadata_MOA":
            moa_meta = _load_opentargets_moa_info()
            adata_for_eval = _merge_with_duplication(adata.copy(), moa_meta)

            # subset adata to drugs with MOA info
            adata_for_eval = adata_for_eval[~adata_for_eval.obs[eval_key].isna()].copy()

        else:
            adata_for_eval = adata.copy()

        if eval_key not in adata_for_eval.obs.columns:
            raise ValueError(f"Eval key '{eval_key}' not in metadata")

        bm = Benchmarker(
            adata_for_eval,
            batch_key=batch_key,
            label_key=eval_key,
            embedding_obsm_keys=methods.split(" "),
        )
        bm.benchmark()

        df = bm.get_results(min_max_scale=False)

        # ugly code to convert to tidy format so we can properly save to parquet
        # TODO(ttreis): maybe refactor this to be more readable
        metric_types = df.loc["Metric Type"].to_frame(name="Metric Type").reset_index()
        metric_types.columns = ["Metric", "Metric_Type"]
        df2 = df.drop("Metric Type").reset_index()
        df2 = df2.rename(columns={"Embedding": "Method"})
        df_long = df2.melt(id_vars=["Method"], var_name="Metric", value_name="Value")
        df_tidy = df_long.merge(metric_types, on="Metric", how="left")
        df_tidy = df_tidy.query("Metric_Type != 'Aggregate score'")
        df_tidy["eval_key"] = eval_key

        results.append(df_tidy)

    # make pretty
    df_tidy = pd.concat(results)
    df_tidy = df_tidy[["eval_key", "Method", "Metric", "Value", "Metric_Type"]]
    df_tidy.columns = df_tidy.columns.str.replace(r"\s+", "_", regex=True).str.lower()
    cols_to_modify = ["method", "metric", "metric_type"]
    df_tidy[cols_to_modify] = df_tidy[cols_to_modify].applymap(
        lambda col: col.replace(" ", "_").lower() if isinstance(col, str) else col
    )

    df_tidy.reset_index(drop=True).to_parquet(output_path, index=False)


if __name__ == "__main__":
    adata_path = sys.argv[1]
    output_path = sys.argv[2]
    batch_key = sys.argv[3]
    eval_keys = sys.argv[4]
    methods = sys.argv[5]
    run_scibmetrics_benchmarker(adata_path, output_path, batch_key, eval_keys, methods)
