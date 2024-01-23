import itertools
from functools import partial

import matplotlib
import pandas as pd
from matplotlib import pyplot as plt
from plottable import ColumnDefinition, Table
from plottable.plots import bar

from plot.colors import (
    BATCH_CMAP,
    COMPOUND_COLORS,
    MICRO_CMAP,
    POSCON_MAP,
    POSCONS,
    SOURCE_CMAP,
)


def despine(ax):
    for side in ["left", "right", "top", "bottom"]:
        ax.spines[side].set_visible(False)
        ax.spines[side].set_visible(False)
    ax.tick_params(which="both",
                   bottom=False,
                   left=False,
                   labelbottom=False,
                   labelleft=False)


def query_multiple_pos(embds: pd.DataFrame):
    multiple_pos = (embds.groupby("Compound")["Metadata_Well"].nunique()
                    [lambda x: x > 1].index)
    multiple_pos = embds[embds["Compound"].isin(multiple_pos)]
    return multiple_pos


def create_compound_cmap(embds: pd.DataFrame):
    compounds = embds["Compound"].drop_duplicates().tolist()
    non_poscons = [c for c in compounds if c not in POSCONS]
    poscons = [c for c in compounds if c in POSCONS]
    order = poscons + non_poscons
    compound_cmap = dict(zip(non_poscons, itertools.cycle(COMPOUND_COLORS)))
    compound_cmap.update(POSCON_MAP)
    return compound_cmap, order


def scatter_panel(embds,
                  fig: plt.Figure,
                  spec: plt.GridSpec,
                  row: int,
                  title=False):
    methods = embds["method"].drop_duplicates().to_list()
    for i, method in enumerate(methods):
        points = embds.query("method==@method").sample(frac=1)
        x, y = points["x"], points["y"]
        colors = points["colors"]
        ax = fig.add_subplot(spec[row, i])
        ax.scatter(x, y, c=colors, s=6)
        despine(ax)
        if title:
            ax.set_title(method)


def white_yellow_green_cm():
    lut_size = 256
    spec = [
        (1.0, 1.0, 1.0),
        (0.90196078431372551, 0.96078431372549022, 0.81568627450980391),
        (0.30196078431372547, 0.5725490196078431, 0.12941176470588237),
    ]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "wYlGn", spec, lut_size)
    return cmap


def get_scalar_mapppable(col_data, norm_type=None):
    if norm_type == "minmax":
        vmin = col_data.min()
        vmax = col_data.max()
    if norm_type == "interquartile":
        # taken from plottable.cmap.normed_cmap
        num_stds = 2.5
        _median, _std = s.median(), s.std()
        vmin = _median - num_stds * _std
        vmax = _median + num_stds * _std
    else:
        vmin, vmax = 0, 1

    cmap = white_yellow_green_cm()
    norm = matplotlib.colors.Normalize(vmin, vmax)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    return m


def results_table(pivot_path: str, ax: plt.Axes, min_max_scale: bool = False):
    """
    Adapted from:
    https://github.com/yoseflab/scib-metrics/blob/0.4.1/src/scib_metrics/benchmark/_core.py#L276-L364
    """
    df = pd.read_parquet(pivot_path)

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
        mappable = get_scalar_mapppable(df.droplevel(0, axis="columns")[col])
        col_def = ColumnDefinition(
            col,
            title=col.replace(" ", "\n", 1),
            textprops=textprops,
            width=1,
            cmap=mappable.to_rgba,
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


def all_panel():
    pivot_path = "outputs/scenario_4/plots/data/pivot_scores.parquet"
    embd_path = "outputs/scenario_4/plots/data/embeddings.parquet"
    embds = pd.read_parquet(embd_path)
    embds = embds.query('~Compound.str.startswith("DMSO")')

    # Sort based on the pivot scores table
    rank = pd.read_parquet(pivot_path).index.to_list()
    # Setting Baseline first
    # rank.remove('Baseline')
    # rank.insert(0 , 'Baseline')

    embds = embds.set_index("method").loc[rank].reset_index()

    fig = plt.figure(figsize=(20, 15))
    spec = fig.add_gridspec(6, 7, height_ratios=[2.5, 0.1, 1, 1, 1, 0.7])

    # Table
    ax = fig.add_subplot(spec[0, :])
    results_table("outputs/scenario_4/plots/data/pivot_scores.parquet", ax)

    cpd_embds = query_multiple_pos(embds).copy()
    compound_cmap, cpd_order = create_compound_cmap(cpd_embds)
    cpd_embds["colors"] = cpd_embds["Compound"].map(compound_cmap)
    scatter_panel(cpd_embds, fig, spec, row=2, title=True)

    embds["colors"] = embds["Source"].map(SOURCE_CMAP)
    scatter_panel(embds, fig, spec, row=3)

    embds["colors"] = embds["Microscope"].map(MICRO_CMAP)
    scatter_panel(embds, fig, spec, row=4)

    ax_cpd = fig.add_subplot(spec[5, 0])
    colors = [compound_cmap[lbl] for lbl in cpd_order]
    hidden_trace = partial(ax_cpd.scatter, x=[], y=[], ls="", marker="o")

    ax_cpd.legend(
        handles=[hidden_trace(color=c) for c in colors],
        labels=cpd_order,
        loc="upper left",
        title="Compound",
        ncols=4,
    )

    ax_micro = fig.add_subplot(spec[5, 4])
    labels = embds["Microscope"].drop_duplicates().to_list()
    colors = [MICRO_CMAP[lbl] for lbl in labels]

    hidden_trace = partial(ax_micro.scatter, x=[], y=[], ls="", marker="o")
    ax_micro.legend(
        handles=[hidden_trace(color=c) for c in colors],
        labels=labels,
        loc="upper left",
        title="Microscope",
        ncols=1,
    )

    ax_source = fig.add_subplot(spec[5, 6])
    labels = embds["Source"].drop_duplicates().to_list()
    colors = [SOURCE_CMAP[lbl] for lbl in labels]

    hidden_trace = partial(ax_source.scatter, x=[], y=[], ls="", marker="o")
    ax_source.legend(
        handles=[hidden_trace(color=c) for c in colors],
        labels=labels,
        loc="upper right",
        title="Source",
        ncols=1,
    )
    despine(ax_source)
    despine(ax_micro)
    despine(ax_cpd)

    ax_colorbar = fig.add_subplot(spec[1, 3:5])
    fig.colorbar(
        get_scalar_mapppable([]),
        cax=ax_colorbar,
        orientation="horizontal",
    )
    ax_colorbar.xaxis.set_ticks_position("top")
    ax_colorbar.set_ylabel("Score ", ha="right", rotation="horizontal")


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


all_panel()
plt.savefig("fig.png", bbox_inches="tight")
