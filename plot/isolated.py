import os
import tempfile

import pandas as pd

from metrics.map import average_precision_negcon, mean_average_precision

model = "seurat_cca"
ap_scores = pd.read_parquet(
    f"./outputs/scenario_3/metrics/prod/mad_int_featselect_{model}_ap_negcon.parquet"
)
map_scores = pd.read_parquet(
    f"./outputs/scenario_3/metrics/prod/mad_int_featselect_{model}_map_negcon.parquet"
)
df = pd.read_parquet(f"./outputs/scenario_3/mad_int_featselect_{model}.parquet")
occ = (
    ap_scores.groupby("Metadata_JCP2022", observed=True)["Metadata_Source"]
    .apply(set)
    .reset_index()
)
sources = ap_scores["Metadata_Source"].unique()
for source in sources:
    occ[source] = occ["Metadata_Source"].apply(lambda x: source in x)
occ = occ.drop(columns="Metadata_Source").set_index("Metadata_JCP2022")
cpds_2_6 = occ.query("source_2 and source_6 and not source_10").index
cpds_2_6_10 = occ.query("source_2 and source_6 and source_10").index

with tempfile.TemporaryDirectory() as tmp:
    ap_path = os.path.join(tmp, "ap.parquet")
    map_path = os.path.join(tmp, "map.parquet")
    parquet_path = os.path.join(tmp, "feats.parquet")
    query = "(Metadata_JCP2022=='DMSO' or Metadata_JCP2022 in @cpds_2_6_10) and Metadata_Source != 'source_10'"
    subset = df.query(query)
    subset.to_parquet(parquet_path, index=False)
    average_precision_negcon(parquet_path, ap_path, plate_types=["COMPOUND"])
    mean_average_precision(ap_path, map_path)
    map_scores_three_sources = pd.read_parquet(map_path)

map_scores_two_sources = map_scores.query("Metadata_JCP2022 in @cpds_2_6")

### PLOT ###
import base64
from io import BytesIO

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def show_inline(close=False):
    iobytes = BytesIO()
    plt.savefig(iobytes, format="png", bbox_inches="tight", dpi=300)
    iobytes.seek(0)
    encode = base64.b64encode(iobytes.read())
    props = "File=inline=1"
    print(f"\x1b]1337;{props}:" + encode.decode() + "\x07\n")
    if close:
        plt.close()


df_plot = pd.DataFrame(
    dict(
        mAP=np.concatenate(
            [
                map_scores_two_sources.mean_average_precision,
                map_scores_three_sources.mean_average_precision,
            ]
        ),
        subset=np.repeat(
            ["two sources", "three sources"],
            [len(map_scores_two_sources), len(map_scores_three_sources)],
        ),
    )
)

fig, axes = plt.subplots(2, sharex=True, gridspec_kw={"height_ratios": (0.15, 0.85)})
sns.boxplot(df_plot, hue="subset", x="mAP", ax=axes[0], gap=0.5)
axes[0].get_legend().set_visible(False)
axes[0].set_yticks([])
axes[0].set_ylabel(None)
sns.histplot(
    df_plot, hue="subset", x="mAP", stat="probability", common_norm=False, ax=axes[1]
)
plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.savefig("figures/sup_figure_F.pdf", format="pdf", bbox_inches="tight", dpi=300)
show_inline(True)
