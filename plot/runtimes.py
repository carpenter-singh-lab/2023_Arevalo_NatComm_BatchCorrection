import base64
import datetime
import json
from glob import glob
from io import BytesIO

import matplotlib
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from metrics import DIMENSION_MAP
from plot.colors import METHOD_FMT, METRIC_FMT
from preprocessing.io import get_num_rows


def load_json(f: str):
    with open(f) as f_in:
        return json.load(f_in)


def timeTicks(x, pos):
    d = datetime.timedelta(seconds=int(x))
    return str(d)


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
    "harmony": "Method",
    "scanorama": "Method",
    "mnn": "Method",
    "sphering_explore": "Method",
    "combat": "Method",
    "INT": "Preprocessing",
    # "compute_negcon_stats": "Preprocessing",
    "select_variant_feats": "Preprocessing",
    "mad_normalize": "Preprocessing",
    # "write_parquet": "Preprocessing",
    "featselect": "Preprocessing",
    "clustering": "bio",
    "average_precision_negcon": "bio",
    "average_precision_nonrep": "bio",
}
categories.update(DIMENSION_MAP)


files = glob(".snakemake/metadata/*")
df = pd.DataFrame(map(load_json, files))
df["runtime"] = df.eval("(endtime - starttime)")
df["category"] = df["rule"].map(categories)
order = df.groupby("rule")["runtime"].mean().sort_values().index
order = pd.Series(index=order, data=range(len(order)))
df = df.sort_values(by="rule", key=order.get).sort_values(by="category", kind="stable")
df["scenario"] = df.input.apply(" ".join).str.extract(r"(scenario_.)", expand=False)
df["rule"] = df.rule.apply(lambda x: METHOD_FMT.get(x, x))
df["rule"] = df.rule.apply(lambda x: METRIC_FMT.get(x, x))
num_comp_map = dict(
    scenario_1="~300",
    scenario_2="~300",
    scenario_3="~80K",
    scenario_4="~300",
    scenario_5="~80K",
)
df["Number of compounds"] = df["scenario"].map(num_comp_map)
num_samp_map = {
    c: get_num_rows(f"./outputs/{c}/raw.parquet") for c in num_comp_map.keys()
}
df["Number of samples"] = df["scenario"].map(num_samp_map)

# Rename categories
df["category"].replace(
    {"batch": "Batch correction", "bio": "Bio metrics"}, inplace=True
)
# Rename rules
df["rule"].replace(
    {
        "sphering_explore": "Sphering",
        "average_precision_nonrep": "mAP",
        "average_precision_negcon": "mAP",
        "select_variant_feats": "Remove low-variance features",
        "INT": "INT transformation",
        "mad_normalize": "MAD normalization",
        "featselect": "Feature selection",
    },
    inplace=True,
)

fig, axes = plt.subplots(
    ncols=4, nrows=1, figsize=(24, 8), constrained_layout=True, sharey=True, sharex=True
)

formatter = matplotlib.ticker.FuncFormatter(timeTicks)
for ax, category in zip(
    axes, ["Method", "Batch correction", "Bio metrics", "Preprocessing"]
):
    subdf = df.query("category==@category")
    if category == "Method":
        subdf = subdf.sort_values(
            by=["scenario", "rule", "runtime"],
            ascending=[True, True, False],
        ).drop_duplicates(["scenario", "rule"])
    g = sns.lineplot(subdf, x="Number of samples", y="runtime", hue="rule", ax=ax)
    g.get_legend().set_title(None)
    ax.set_title(category)
    ax.set_ylabel("Runtime (hh:mm:ss)")
    ax.set_xlabel(None)
    ax.grid(axis="y", zorder=0)
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.yaxis.set_major_formatter(formatter)
    ax.yaxis.set_ticks([], minor=True)
    ax.yaxis.set_ticks([10, 60, 600, 3600, 3600 * 2, 3600 * 5, 3600 * 10])
    ax.xaxis.set_ticks([], minor=True)
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_ticks(list(num_samp_map.values()))
    ax.set_xticklabels(list(num_samp_map.values()), rotation=40, ha="right")
fig.supxlabel("Number of samples")
show_inline(True)
