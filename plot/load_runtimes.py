import base64
import json
from glob import glob
from io import BytesIO

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def show_inline(close=False):
    iobytes = BytesIO()
    plt.savefig(iobytes, format="png", bbox_inches="tight", dpi=300)
    iobytes.seek(0)
    encode = base64.b64encode(iobytes.read())
    props = "File=inline=1"
    print(f"\x1b]1337;{props}:" + encode.decode() + "\x07\n")
    if close:
        plt.close()


categories = {
    "nmi": "metric",
    "ari": "metric",
    "mean_average_precision": "metric",
    "pcr": "metric",
    "average_precision_negcon": "metric",
    "compute_negcon_stats": "preprocessing",
    "combat": "method",
    "select_variant_feats": "preprocessing",
    "pcr_batch": "metric",
    "average_precision_nonrep": "metric",
    "lisi_batch": "metric",
    "sphering_explore": "method",
    "mad_normalize": "preprocessing",
    "write_parquet": "preprocessing",
    "silhouette_batch": "metric",
    "featselect": "preprocessing",
    "INT": "preprocessing",
    "clustering": "metric",
    "graph_conn": "metric",
    "asw": "metric",
    "harmony": "method",
    "lisi_label": "metric",
    "mnn": "method",
}


def load_json(f: str):
    with open(f) as f_in:
        return json.load(f_in)


files = glob(".snakemake/metadata/*")
df = pd.DataFrame(map(load_json, files))
df["runtime"] = df.eval("endtime - starttime")
df["category"] = df["rule"].map(categories)
order = df.groupby("rule")["runtime"].mean().sort_values().index
order = pd.Series(index=order, data=range(len(order)))
df = (
    df.sort_values(by="rule", key=order.get)
    .sort_values(by="category", kind="stable")
    .dropna(subset="category")
)
g = sns.barplot(df, x="rule", y="runtime", hue="category")
g.set_yscale("log")
plt.xticks(rotation=45, ha="right")
show_inline(True)
