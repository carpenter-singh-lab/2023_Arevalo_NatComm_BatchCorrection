import warnings
from difflib import SequenceMatcher

import numpy as np
import pandas as pd
from sklearn.preprocessing import minmax_scale

from .colors import METHOD_FMT, METRIC_FMT

# METHOD_FMT = {
#     'combat': 'Combat',
#     'harmony': 'Harmony',
#     'scvi': 'scVI',
#     'sysvi': 'sysVI',
#     'scpoli': 'scPoli',
#     'scanorama': 'Scanorama',
#     'mnn': 'MNN',
#     'desc': 'DESC',
#     'baseline': 'Baseline',
#     'sphering': 'Sphering',
#     'fastMNN': 'fastMNN',
#     'seurat_rpca': 'Seurat RPCA',
#     'seurat_cca': 'Seurat CCA'
# }

# METRIC_FMT = {
#     'silhouette_batch': 'Silhouette batch',
#     'pcr_batch': 'PCR comparison',
#     'pcr': 'PCR',
#     'graph_conn': 'Graph connectivity',
#     'kbet': 'KBET',
#     'lisi_batch': 'LISI batch',
#     'lisi_label': 'LISI label',
#     'negcon_mean_map': 'mAP (controls)',
#     'negcon_fraction_below_p': '%Retrieved (controls)',
#     'negcon_fraction_below_corrected_p': '%Retrieved-corr (controls)',
#     'nonrep_mean_map': 'mAP (nonrep)',
#     'nonrep_fraction_below_p': '%Retrieved (nonrep)',
#     'nonrep_fraction_below_corrected_p': '%Retrieved-corr (nonrep)',
#     'nmi': 'Leiden NMI',
#     'ari': 'Leiden ARI',
#     'asw': 'Silhouette label',
#     'il_f1': 'Isolated labels(F1)',
#     'il_asw': 'Isolated labels(ASW)',
# }

def _common_prefix_suffix(strings: list[str]):
    prefix = strings[0]
    suffix = prefix[::-1]
    for string in strings[1:]:
        match = SequenceMatcher(a=prefix, b=string).get_matching_blocks()[0]
        prefix = prefix[:match.size]
        match = SequenceMatcher(a=suffix,
                                b=string[::-1]).get_matching_blocks()[0]
        suffix = suffix[:match.size]
    suffix = suffix[::-1]
    return prefix, suffix


def _jitter(trace, win_size=0.02):
    range_x = max(trace.x) - min(trace.x)
    factor = range_x * win_size
    trace.x += np.random.uniform(-factor, factor, [len(trace.x)])

    range_y = max(trace.y) - min(trace.y)
    factor = range_y * win_size
    trace.y += np.random.uniform(-factor, factor, [len(trace.y)])


def load_all_parquet(files, key_name="file_id", placeholder="baseline"):
    dframe = []
    dframes = [pd.read_parquet(file) for file in files]
    prefix, suffix = _common_prefix_suffix(files)
    start = len(prefix)
    end = -len(suffix)
    for dframe, file in zip(dframes, files):
        dframe[key_name] = file[start:end]
    dframe = pd.concat(dframes)
    dframe[key_name] = dframe[key_name].str.strip("_").replace("", placeholder)
    return dframe


def prepare_embeddings(embd_files: list[str], output_path: str, anon=True):
    embds = load_all_parquet(embd_files, key_name="method")
    embds["method"] = embds["method"].map(lambda x: METHOD_FMT.get(x, x))

    # jitter embds
    embds_jitter = []
    for _, embd in embds.groupby("method"):
        _jitter(embd)
        with warnings.catch_warnings():
            warnings.simplefilter(action="ignore", category=FutureWarning)
            # hack to make plotly happy in ranges
            embd[["x", "y"]] = minmax_scale(embd[["x", "y"]])
        embds_jitter.append(embd)
    embds = pd.concat(embds_jitter)

    if anon:
        # Make batches and sources annonymous
        renamer = {
            batch: f"{i:02d}"
            for i, batch in enumerate(embds["Metadata_Batch"].unique(), 1)
        }
        embds["Metadata_Batch"] = embds["Metadata_Batch"].apply(renamer.get)

        renamer = {
            source: f'{int(source.split("_")[-1]):02d}'
            for source in embds["Metadata_Source"].unique()
        }
        embds["Metadata_Source"] = embds["Metadata_Source"].apply(renamer.get)

    embds = embds.rename(
        columns={
            "Metadata_Batch": "Batch",
            "Metadata_Microscope": "Microscope",
            "Metadata_Source": "Source",
            "Metadata_JCP2022": "Compound",
            "Metadata_Row": "Row",
            "Metadata_Column": "Column",
            "method": "Method",
        })
    embds.to_parquet(output_path, index=False)


def query_multiple_pos(embds: pd.DataFrame):
    multiple_pos = (embds.groupby("Compound")["Metadata_Well"].nunique()
                    [lambda x: x > 1].index)
    multiple_pos = embds[embds["Compound"].isin(multiple_pos)]
    return multiple_pos


def tidy_scores(metrics_files, metrics_redlist, methods_redlist, tidy_path):
    scores = load_all_parquet(metrics_files, key_name="method")
    scores = scores.query("metric not in @metrics_redlist")
    scores = scores[scores["method"].apply(
        lambda x: all(m not in x for m in methods_redlist))]
    scores.to_parquet(tidy_path, index=False)


def pivot_scores(tidy_path, pivot_path):
    scores = pd.read_parquet(tidy_path)
    scores = scores.pivot_table(
        index="method",
        columns=["eval_key", "metric_type", "metric"],
        values="value"
    )
    scores.index.name = "Method"
    scores.to_parquet(pivot_path)
