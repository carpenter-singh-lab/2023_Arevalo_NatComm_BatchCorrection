import pandas as pd
from copairs.map import mean_average_precision
from copairs.map.multilabel import average_precision

from preprocessing import io
from preprocessing.metadata import MAPPER


def get_labels(drugs_url, inchi_mapper_url, compound_path, output_path):
    drugs = pd.read_csv(drugs_url, skiprows=9, sep="\t").dropna(subset="target")
    mapper = pd.read_csv(inchi_mapper_url, sep="\t").drop(columns="target")
    drugs["pert_iname"] = drugs["pert_iname"].str.lower()
    mapper["pert_iname"] = mapper["pert_iname"].str.lower()
    targets = mapper.merge(drugs, on="pert_iname")
    targets = targets.groupby("InChIKey")["target"].apply("|".join)
    targets = targets.apply(lambda x: list(set(x.split("|"))))

    compound = pd.read_csv(compound_path)
    compound["Metadata_JCP2022"] = compound["Metadata_JCP2022"].apply(
        lambda x: MAPPER.get(x, x)
    )
    compound["Metadata_Target"] = compound["Metadata_InChIKey"].map(targets)

    labels = compound[["Metadata_JCP2022", "Metadata_Target"]].dropna()
    labels.to_parquet(output_path, index=False)


def median_profile(dframe_path, median_path):
    dframe = pd.read_parquet(dframe_path)
    featcols = [c for c in dframe.columns if not c.startswith("Meta")]
    agg_funcs = {c: "median" for c in featcols}
    dframe = dframe.groupby("Metadata_JCP2022", observed=True).agg(agg_funcs)
    dframe.reset_index().to_parquet(median_path)


def annotate_median_profile(dframe_path, annotations_path, output_path):
    annotations = pd.read_parquet(annotations_path)
    dframe = pd.read_parquet(dframe_path)
    mapper = annotations.set_index("Metadata_JCP2022")["Metadata_Target"]
    dframe["Metadata_Target"] = dframe["Metadata_JCP2022"].map(mapper)
    dframe = dframe.dropna(subset="Metadata_Target")
    dframe.to_parquet(output_path, index=False)


def target_ap(dframe_path, ap_path):
    meta, vals, _ = io.split_parquet(dframe_path)
    result = average_precision(
        meta,
        vals,
        pos_sameby=["Metadata_Target"],
        pos_diffby=["Metadata_JCP2022"],
        neg_sameby=[],
        neg_diffby=["Metadata_Target"],
        batch_size=20000,
        multilabel_col="Metadata_Target",
    )
    result.reset_index(drop=True).to_parquet(ap_path)


def target_map(ap_path, map_path, threshold=0.05):
    ap_scores = pd.read_parquet(ap_path)
    map_scores = mean_average_precision(
        ap_scores, "Metadata_Target", threshold=threshold, null_size=10000, seed=0
    )
    map_scores.to_parquet(map_path)
