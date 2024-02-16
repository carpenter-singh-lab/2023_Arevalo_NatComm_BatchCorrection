import base64
import datetime
import json
from glob import glob
from io import BytesIO

import matplotlib
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from plot.colors import METHOD_FMT, METRIC_FMT


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
    "nmi": "Metric",
    "ari": "Metric",
    "mean_average_precision": "Metric",
    "pcr": "Metric",
    "average_precision_negcon": "Metric",
    "compute_negcon_stats": "Preprocessing",
    "combat": "Method",
    "select_variant_feats": "Preprocessing",
    "pcr_batch": "Metric",
    "average_precision_nonrep": "Metric",
    "lisi_batch": "Metric",
    "sphering_explore": "Method",
    "mad_normalize": "Preprocessing",
    "write_parquet": "Preprocessing",
    "silhouette_batch": "Metric",
    "featselect": "Preprocessing",
    "INT": "Preprocessing",
    "clustering": "Metric",
    "graph_conn": "Metric",
    "asw": "Metric",
    "harmony": "Method",
    "lisi_label": "Metric",
    "mnn": "Method",
    "kbet": "Metric",
    "select_best_sphering": "",
    "umap": "Visualization",
}


def load_json(f: str):
    with open(f) as f_in:
        return json.load(f_in)


files = glob(".snakemake/metadata/*")
df = pd.DataFrame(map(load_json, files))
df["runtime"] = df.eval("(endtime - starttime)")
df["category"] = df["rule"].map(categories)
order = df.groupby("rule")["runtime"].mean().sort_values().index
order = pd.Series(index=order, data=range(len(order)))
df = (
    df.sort_values(by="rule", key=order.get).sort_values(by="category", kind="stable")
    # .dropna(subset="category")
)
df["Scenario"] = df.input.apply(" ".join).str.extract(r"(scenario_.)", expand=False)
df["rule"] = df.rule.apply(lambda x: METHOD_FMT.get(x, x))
df["rule"] = df.rule.apply(lambda x: METRIC_FMT.get(x, x))
ax = sns.barplot(df, x="rule", y="runtime", hue="category", zorder=2)
ax.grid(axis="y", zorder=0)
ax.set_yscale("log")


def timeTicks(x, pos):
    d = datetime.timedelta(seconds=int(x))
    return str(d)


ax.yaxis.set_ticks([10, 60, 600, 3600, 3600 * 2, 3600 * 5])
formatter = matplotlib.ticker.FuncFormatter(timeTicks)
ax.yaxis.set_major_formatter(formatter)
ax.yaxis.set_ticks([], minor=True)
ax.set_ylabel("Runtime (hh:mm:ss)")
ax.set_xlabel("Operation")
ax.get_legend().set_title("Category")

plt.xticks(rotation=45, ha="right")
show_inline(True)
