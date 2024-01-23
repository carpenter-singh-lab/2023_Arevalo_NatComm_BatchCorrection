import itertools
from functools import partial

import pandas as pd
from matplotlib import pyplot as plt
from plot.table import draw as draw_table, add_colorbar

from plot.colors import (
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


def all_panel(embd_path: str, pivot_path: str):
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

    # Table with colorbar
    ax = fig.add_subplot(spec[0, :])
    draw_table(pivot_path, ax)
    ax_colorbar = fig.add_subplot(spec[1, 3:5])
    add_colorbar(fig, ax_colorbar)

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



pivot_path = "outputs/scenario_4/plots/data/pivot_scores.parquet"
embd_path = "outputs/scenario_4/plots/data/embeddings.parquet"
all_panel(embd_path, pivot_path)
plt.savefig("fig.png", bbox_inches="tight")
