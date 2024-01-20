import itertools

import matplotlib
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from plottable import ColumnDefinition, Table
from plottable.cmap import normed_cmap
from plottable.plots import bar

from plot.colors import (
    BATCH_CMAP,
    COMPOUND_COLORS,
    METHOD_FMT,
    METRIC_FMT,
    MICRO_CMAP,
    POSCON_MAP,
    POSCONS,
    SOURCE_CMAP,
)


def cmap_table_fn(col_data):
    return normed_cmap(col_data, cmap=matplotlib.cm.PRGn, num_stds=2.5)


def results_table(pivot_path: str, ax, min_max_scale: bool = False):
    """
    Adapted from:
    https://github.com/yoseflab/scib-metrics/blob/0.4.1/src/scib_metrics/benchmark/_core.py#L276-L364
    """
    df = pd.read_parquet(pivot_path)
    if min_max_scale:
        df = pd.DataFrame(minmax_scale(df), index=df.index, columns=df.columns)
        cols = ["Batch correction", "Bio metrics"]
        for col in cols:
            df["mean", col] = df[col].mean(axis=1)
        total = np.array([0.4, 0.6]) @ df["mean"][cols].T
        df["mean", "Total"] = total
        df = df.sort_values(by=("mean", "Total"), ascending=False)

    column_definitions = [
        ColumnDefinition("method",
                         width=1.5,
                         textprops={
                             "ha": "left",
                             "weight": "bold"
                         }),
    ]
    score_cols = df[["Batch correction",
                     "Bio metrics"]].columns.get_level_values(1)
    textprops = {"ha": "center", "bbox": {"boxstyle": "circle", "pad": 0.25}}
    groupmap = dict(df.columns.swaplevel())
    for i, col in enumerate(score_cols):
        col_def = ColumnDefinition(
            col,
            title=col.replace(" ", "\n", 1),
            textprops=textprops,
            width=1,
            cmap=cmap_table_fn(df.droplevel(0, axis="columns")[col]),
            group=groupmap[col],
            formatter="{:.2f}",
        )
        column_definitions.append(col_def)

    plot_kw = {
        "cmap": matplotlib.cm.YlGnBu,
        "plot_bg_bar": False,
        "annotate": True,
        "height": 0.9,
        "formatter": "{:.2f}",
    }
    agg_cols = df["mean"].columns
    for i, col in enumerate(agg_cols):
        col_def = ColumnDefinition(
            col,
            width=1,
            plot_kw=plot_kw,
            title=col.replace(" ", "\n", 1),
            plot_fn=bar,
            group="Aggregate score",
            border="left" if i == 0 else None,
        )
        column_definitions.append(col_def)

    plt.style.use("default")
    # fig, ax = plt.subplots(figsize=(len(df.columns) * 1.25, 3 + 0.3 * len(df)))
    tab = Table(
        df.droplevel(0, axis="columns").reset_index(),
        cell_kw={
            "linewidth": 0,
            "edgecolor": "k",
        },
        column_definitions=column_definitions,
        ax=ax,
        row_dividers=True,
        footer_divider=True,
        textprops={
            "fontsize": 10,
            "ha": "center"
        },
        row_divider_kw={
            "linewidth": 1,
            "linestyle": (0, (1, 5))
        },
        col_label_divider_kw={
            "linewidth": 1,
            "linestyle": "-"
        },
        column_border_kw={
            "linewidth": 1,
            "linestyle": "-"
        },
        index_col="method",
    )
    tab.autoset_fontcolors(colnames=list(df.columns.get_level_values(1)))


def annotate_axes(ax, text, fontsize=18):
    ax.text(
        0.5,
        0.5,
        text,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=fontsize,
        color="darkgrey",
    )


embds = pd.read_parquet("outputs/scenario_1/plots/data/embeddings.parquet")
embds = embds.query('~Compound.str.startswith("DMSO")')

fig = plt.figure(figsize=(20, 13))
spec = fig.add_gridspec(4, 7, height_ratios=[2.3, 1, 1, 0.7])

ax = fig.add_subplot(spec[0, :])
results_table("outputs/scenario_1/plots/data/pivot_scores.parquet", ax)
method_order = embds["method"].drop_duplicates()


def despine(ax):
    for side in ["left", "right", "top", "bottom"]:
        ax.spines[side].set_visible(False)
        ax.spines[side].set_visible(False)
    ax.tick_params(which="both",
                   bottom=False,
                   left=False,
                   labelbottom=False,
                   labelleft=False)


def plot_umap(ax, x, y, colors):
    ax.scatter(x, y, c=colors, s=6)
    despine(ax)


multiple_pos = (embds.groupby("Compound")["Metadata_Well"].nunique()
                [lambda x: x > 1].index)
multiple_pos = embds[embds["Compound"].isin(multiple_pos)]
compounds = multiple_pos["Compound"].drop_duplicates().tolist()
non_poscons = [c for c in compounds if c not in POSCONS]
poscons = [c for c in compounds if c in POSCONS]
hue_order = poscons + non_poscons
compound_cmap = dict(zip(non_poscons, itertools.cycle(COMPOUND_COLORS)))
compound_cmap.update(POSCON_MAP)

for i, method in enumerate(method_order):
    group = embds.query("method==@method").sample(frac=1)
    x, y = group["x"], group["y"]
    colors = group["Batch"].map(BATCH_CMAP)
    ax_batch = fig.add_subplot(spec[2, i])
    plot_umap(ax_batch, x, y, colors)

    group = multiple_pos.query("method==@method").sample(frac=1)
    x, y = group["x"], group["y"]
    colors = group["Compound"].map(compound_cmap)
    ax_bio = fig.add_subplot(spec[1, i])
    plot_umap(ax_bio, x, y, colors)

ax_bio = fig.add_subplot(spec[3, 0])
labels = hue_order
colors = [compound_cmap[lbl] for lbl in labels]


def hidden_trace(c):
    return ax_bio.scatter([], [], color=c, ls="", marker="o")


ax_bio.legend(handles=[hidden_trace(c) for c in colors],
              labels=labels,
              loc="upper left",
              ncols=4)

ax_batch = fig.add_subplot(spec[3, 6])
labels = embds["Batch"].drop_duplicates().to_list()
colors = [BATCH_CMAP[lbl] for lbl in labels]


def hidden_trace(c):
    return ax_batch.scatter([], [], color=c, ls="", marker="o")


ax_batch.legend(handles=[hidden_trace(c) for c in colors],
                labels=labels,
                loc="upper right",
                ncols=4)
despine(ax_batch)
despine(ax_bio)
plt.savefig("fig.pdf", bbox_inches="tight")
