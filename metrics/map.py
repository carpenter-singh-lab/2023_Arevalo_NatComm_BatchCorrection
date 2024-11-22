import logging

from typing import Literal
import copairs.map as copairs
import pandas as pd
import anndata as ad

from preprocessing.io import split_parquet
from metrics.scib import _load_opentargets_moa_info, _merge_with_duplication

logger = logging.getLogger(__name__)

def _get_idx_of_samples_for_computation(meta, plate_type, ignore_dmso=False):
    index = meta["Metadata_PlateType"].isin(plate_type)
    index &= (meta["Metadata_PertType"] != "poscon")
    valid_cmpd = meta.loc[index, "Metadata_JCP2022"].value_counts()
    valid_cmpd = valid_cmpd[valid_cmpd > 1].index
    index &= meta["Metadata_JCP2022"].isin(valid_cmpd)
    # TODO: This compound has many more replicates than any other. ignoring it
    # for now. This filter should be done early on.
    index &= (meta["Metadata_JCP2022"] != "JCP2022_033954")
    if ignore_dmso:
        index &= (meta["Metadata_JCP2022"] != "DMSO")
    return meta.index[index]


def _group_negcons(meta: pd.DataFrame):
    """
    Hack to avoid mAP computation for negcons. Assign a unique id for every
    negcon so that no pairs are found for such samples.
    """
    negcon_ix = (meta["Metadata_JCP2022"] == "DMSO")
    n_negcon = negcon_ix.sum()
    negcon_ids = [f"DMSO_{i}" for i in range(n_negcon)]
    pert_id = meta["Metadata_JCP2022"].astype("category").cat.add_categories(negcon_ids)
    pert_id[negcon_ix] = negcon_ids
    meta["Metadata_JCP2022"] = pert_id


# def average_precision_negcon(parquet_path, ap_path, plate_types):
#     meta, vals, _ = split_parquet(parquet_path)
#     ix = _get_idx_of_samples_for_computation(meta, plate_types)
#     meta = meta[ix].copy()
#     vals = vals[ix]
#     _group_negcons(meta)
#     result = copairs.average_precision(
#         meta,
#         vals,
#         pos_sameby=["Metadata_JCP2022"],
#         pos_diffby=["Metadata_Well"],
#         neg_sameby=["Metadata_Plate"],
#         neg_diffby=["Metadata_PertType", "Metadata_JCP2022"],
#         batch_size=20000)
#     result = result.query("Metadata_PertType!="negcon"")
#     result.reset_index(drop=True).to_parquet(ap_path)


# def average_precision_nonrep(parquet_path, ap_path, plate_types):
#     meta, vals, _ = split_parquet(parquet_path)
#     ix = _get_idx_of_samples_for_computation(meta, plate_types, ignore_dmso=True)
#     meta = meta[ix].copy()
#     vals = vals[ix]
#     result = copairs.average_precision(
#         meta,
#         vals,
#         pos_sameby=["Metadata_JCP2022"],
#         pos_diffby=[],
#         neg_sameby=["Metadata_Plate"],
#         neg_diffby=["Metadata_JCP2022"],
#         batch_size=20000,
#     )
#     result.reset_index(drop=True).to_parquet(ap_path)

def calculate_map_summary(map_scores: pd.DataFrame):
    map_scores = map_scores.dropna()
    frac_p = map_scores["below_p"].sum() / len(map_scores)
    frac_q = map_scores["below_corrected_p"].sum() / len(map_scores)
    mean_map = map_scores["mean_average_precision"].mean()
    return mean_map, frac_p, frac_q

def calculate_map(
    mode: Literal["negcon", "nonrep"],
    metadata: pd.DataFrame,
    embedding: pd.DataFrame,
    plate_type: str,
    eval_key: str,
    threshold=0.05,
):
    ignore_dmso = mode == "nonrep"
    idx = _get_idx_of_samples_for_computation(
        meta=metadata,
        plate_type=plate_type, 
        ignore_dmso=ignore_dmso
    )
    
    meta = metadata.loc[idx]
    vals = embedding.loc[idx].values # copairs needs np.array
    
    if mode == "negcon":
        _group_negcons(meta)
    
    ap_scores = copairs.average_precision(
        meta,
        vals,
        pos_sameby=[eval_key],
        pos_diffby=[] if mode == "nonrep" else ["Metadata_Well"],
        neg_sameby=["Metadata_Plate"],
        neg_diffby=[eval_key] if mode == "nonrep" else ["Metadata_PertType", eval_key],
        batch_size=20000
    )
    
    if mode == "negcon":
        ap_scores = ap_scores.query("Metadata_PertType!='negcon'")
    
    map_scores = copairs.mean_average_precision(
        ap_scores,
        eval_key,
        threshold=threshold,
        null_size=10000,
        seed=0
    )
    
    mean_map, frac_p, frac_q = calculate_map_summary(map_scores)

    return {
        "mean_map": mean_map,
        "frac_p": frac_p,
        "frac_q": frac_q,
    }


def mean_average_precision(
    adata_path: str, 
    output_path: str, 
    plate_type: str, 
    eval_keys: list[str],
    threshold=0.05
):
    adata = ad.read_h5ad(adata_path)
    results = {}

    # filter out eval_keys for which we cannot perform mAP calculations
    eval_keys.remove("Metadata_MOA")

    # iterate over all keys x embeddings and calculate mAP scores
    for eval_key in eval_keys:
        results[eval_key] = {}
        for embedding in adata.obsm.keys():
            results[eval_key][embedding.lower()] = {}
            for mode in ["negcon", "nonrep"]:
            
                logger.info(f"Calculating {mode} mAP with eval_key '{eval_key}' for '{embedding}'.")
            
                results[eval_key][embedding.lower()][mode] = calculate_map(
                    mode=mode,
                    metadata=adata.obs.copy(),
                    embedding=adata.obsm[embedding].copy(),
                    plate_type=plate_type,
                    eval_key=eval_key,
                    threshold=threshold,
                )
    
    # reshape into tidy format
    rows = []
    for eval_key, methods in results.items():
        for method, modes in methods.items():
            for mode, metrics in modes.items():
                for metric, value in metrics.items():
                    rows.append({
                        "eval_key": eval_key,
                        "method": method,
                        "metric": f"{mode}_{metric}",
                        "value": value,
                    })

    tidy_results = pd.DataFrame(rows)
    tidy_results["metric_type"] = "bio_conservation"
       
    tidy_results.reset_index(drop=True).to_parquet(output_path, index=False)

    logger.info(f"Finished calculating mAP scores for all {len(methods.keys())} embeddings.")