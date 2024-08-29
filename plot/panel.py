import itertools
from functools import partial

import pandas as pd
from matplotlib import pyplot as plt

from .colors import (
    BATCH_CMAP,
    COMPOUND_COLORS,
    MICRO_CMAP,
    POSCON_MAP,
    POSCONS,
    SOURCE_CMAP,
)
from .data import query_multiple_pos
from .scatter import despine, scatter_panel
from .table import add_colorbar
from .table import draw as draw_table


def load_embeddings(embd_path: str, pivot_path: str):
    """Load embeddings filtering DMSO and sorting them by the pivot rank"""
    embds = pd.read_parquet(embd_path)
    embds = embds.query('~Compound.str.startswith("DMSO")')
    rank = pd.read_parquet(pivot_path).index.to_list()
    embds = embds.set_index("Method").loc[rank].reset_index()
    return embds


def add_number(ax: plt.Axes, number: str, **kwargs):
    """Add number in top-left"""
    text_params = dict(
        x=0,
        y=1,
        fontsize=16,
        transform=ax.transAxes,
        fontweight="bold",
        va="center",
        ha="center",
    )
    text_params.update(kwargs)
    ax.text(s=number, **text_params)


def add_table(pivot_path, fig, spec) -> plt.Axes:
    # Table with colorbar
    ax = fig.add_subplot(spec[0, :])
    draw_table(pivot_path, ax)
    ax_colorbar = fig.add_subplot(spec[1, 4:7])
    add_colorbar(fig, ax_colorbar)
    return ax


def add_legend(ax: plt.Axes, cmap: pd.Series, title: str):
    # Remove JCP2022_ prefix
    cmap.index = cmap.index.str.replace("JCP2022_", "")
    hidden_trace = partial(ax.scatter, x=[], y=[], ls="", marker="o")
    # ceiling_division
    ncols = -(len(cmap) // -5)
    ax.legend(
        handles=[hidden_trace(color=c) for c in cmap.values],
        labels=cmap.index.to_list(),
        loc="upper left",
        title=title,
        ncols=ncols,
    )
    despine(ax)


def compound_cmap(embds: pd.DataFrame):
    compounds = query_multiple_pos(embds)["Compound"].drop_duplicates()
    non_poscons = sorted([c for c in compounds if c not in POSCONS])
    poscons = [c for c in compounds if c in POSCONS]
    order = poscons + non_poscons
    cmap = dict(zip(non_poscons, itertools.cycle(COMPOUND_COLORS)))
    cmap.update(POSCON_MAP)
    return order, cmap


def source_cmap(embds: pd.DataFrame):
    order = embds["Source"].drop_duplicates().astype(str).sort_values().to_list()
    return order, SOURCE_CMAP


def micro_cmap(embds: pd.DataFrame):
    order = embds["Microscope"].drop_duplicates().sort_values().to_list()
    return order, MICRO_CMAP


def batch_cmap(embds: pd.DataFrame):
    order = embds["Batch"].drop_duplicates().astype(str).sort_values().to_list()
    return order, BATCH_CMAP


def colorby(embds: pd.DataFrame, column: str) -> pd.Series:
    if column == "Compound":
        order, cmap = compound_cmap(embds)
    elif column == "Source":
        order, cmap = source_cmap(embds)
    elif column == "Microscope":
        order, cmap = micro_cmap(embds)
    elif column == "Batch":
        order, cmap = batch_cmap(embds)

    embds["colors"] = embds[column].map(cmap)
    return pd.Series(cmap)[order]


def results_table(pivot_path: str, fig_path: str):
    fig = plt.figure(figsize=(14, 9))
    spec = fig.add_gridspec(2, 11, height_ratios=[2.5, 0.1])
    add_table(pivot_path, fig, spec)
    plt.savefig(fig_path, bbox_inches="tight")


def full_panel(embd_path: str, pivot_path: str, fig_path: str, scenario: str):
    panel_fn = globals()[scenario]
    panel_fn(embd_path, pivot_path, fig_path)


def scenario_1(embd_path: str, pivot_path: str, fig_path: str):
    fig = plt.figure(figsize=(24, 18))
    spec = fig.add_gridspec(5, 11, height_ratios=[2.5, 0.1, 0.8, 0.8, 0.7])
    ax = add_table(pivot_path, fig, spec)
    add_number(ax, "A", y=0.95)

    embds = load_embeddings(embd_path, pivot_path)

    cmap = colorby(embds, "Compound")
    add_legend(fig.add_subplot(spec[4, 2]), cmap, "Compound")
    axs = scatter_panel(embds.dropna(subset="colors"), fig, spec, row=2, title=True)
    axs[0].set_ylabel("Compound", fontsize=16)
    add_number(axs[0], "B")

    cmap = colorby(embds, "Batch")
    add_legend(fig.add_subplot(spec[4, 7]), cmap, "Batch")
    axs = scatter_panel(embds.dropna(subset="colors"), fig, spec, row=3)
    axs[0].set_ylabel("Batch", fontsize=16)
    add_number(axs[0], "C")

    plt.savefig(fig_path, bbox_inches="tight", dpi=300)


def scenario_2(embd_path: str, pivot_path: str, fig_path: str):
    fig = plt.figure(figsize=(24, 18))
    spec = fig.add_gridspec(5, 11, height_ratios=[2.5, 0.1, 0.8, 0.8, 0.7])
    ax = add_table(pivot_path, fig, spec)
    add_number(ax, "A", y=0.95)

    embds = load_embeddings(embd_path, pivot_path)

    cmap = colorby(embds, "Compound")
    add_legend(fig.add_subplot(spec[4, 2]), cmap, "Compound")
    axs = scatter_panel(embds.dropna(subset="colors"), fig, spec, row=2, title=True)
    axs[0].set_ylabel("Compound", fontsize=16)
    add_number(axs[0], "B")

    cmap = colorby(embds, "Source")
    add_legend(fig.add_subplot(spec[4, 7]), cmap, "Source")
    axs = scatter_panel(embds.dropna(subset="colors"), fig, spec, row=3)
    axs[0].set_ylabel("Source", fontsize=16)
    add_number(axs[0], "C")

    plt.savefig(fig_path, bbox_inches="tight", dpi=300)


def scenario_3(embd_path: str, pivot_path: str, fig_path: str):
    fig = plt.figure(figsize=(24, 18))
    spec = fig.add_gridspec(4, 11, height_ratios=[2, 0.1, 0.6, 0.4])
    ax = add_table(pivot_path, fig, spec)
    add_number(ax, "A", y=0.95)

    embds = load_embeddings(embd_path, pivot_path)

    cmap = colorby(embds, "Source")
    add_legend(fig.add_subplot(spec[3, 5]), cmap, "Source")
    axs = scatter_panel(embds.dropna(subset="colors"), fig, spec, row=2, title=True)
    axs[0].set_ylabel("Source", fontsize=16)
    add_number(axs[0], "B")

    plt.savefig(fig_path, bbox_inches="tight", dpi=300)


def scenario_4(embd_path: str, pivot_path: str, fig_path: str):
    fig = plt.figure(figsize=(24, 18))
    spec = fig.add_gridspec(6, 11, height_ratios=[2.9, 0.1, 1, 1, 1, 0.7])
    ax = add_table(pivot_path, fig, spec)
    add_number(ax, "A", y=0.95)

    embds = load_embeddings(embd_path, pivot_path)

    cmap = colorby(embds, "Compound")
    add_legend(fig.add_subplot(spec[5, 2]), cmap, "Compound")
    axs = scatter_panel(embds.dropna(subset="colors"), fig, spec, row=2, title=True)
    axs[0].set_ylabel("Compound", fontsize=16)
    add_number(axs[0], "B")

    cmap = colorby(embds, "Source")
    add_legend(fig.add_subplot(spec[5, 8]), cmap, "Source")
    axs = scatter_panel(embds.dropna(subset="colors"), fig, spec, row=3)
    axs[0].set_ylabel("Source", fontsize=16)
    add_number(axs[0], "C")

    cmap = colorby(embds, "Microscope")
    add_legend(fig.add_subplot(spec[5, 6]), cmap, "Microscope")
    axs = scatter_panel(embds, fig, spec, row=4)
    axs[0].set_ylabel("Microscope", fontsize=16)
    add_number(axs[0], "D")

    plt.savefig(fig_path, bbox_inches="tight", dpi=300)


def scenario_5(embd_path: str, pivot_path: str, fig_path: str):
    fig = plt.figure(figsize=(24, 18))
    spec = fig.add_gridspec(5, 11, height_ratios=[2.5, 0.1, 0.8, 0.8, 0.7])
    ax = add_table(pivot_path, fig, spec)
    add_number(ax, "A", y=0.95)

    embds = load_embeddings(embd_path, pivot_path)

    cmap = colorby(embds, "Source")
    add_legend(fig.add_subplot(spec[4, 6]), cmap, "Source")
    axs = scatter_panel(embds.dropna(subset="colors"), fig, spec, row=2, title=True)
    axs[0].set_ylabel("Source", fontsize=16)
    add_number(axs[0], "B")

    cmap = colorby(embds, "Microscope")
    add_legend(fig.add_subplot(spec[4, 4]), cmap, "Microscope")
    axs = scatter_panel(embds, fig, spec, row=3)
    axs[0].set_ylabel("Microscope", fontsize=16)
    add_number(axs[0], "C")

    plt.savefig(fig_path, bbox_inches="tight", dpi=300)

def scenario_6(embd_path: str, pivot_path: str, fig_path: str):
    fig = plt.figure(figsize=(24, 18))
    spec = fig.add_gridspec(
        5, 
        13, # adjust when more methods are added
        height_ratios=[2.5, 0.1, 0.8, 0.8, 0.7]
    )
    ax = add_table(pivot_path, fig, spec)
    add_number(ax, "A", y=0.95)

    embds = load_embeddings(embd_path, pivot_path)

    cmap = colorby(embds, "Compound")
    add_legend(fig.add_subplot(spec[4, 2]), cmap, "Compound")
    axs = scatter_panel(embds.dropna(subset="colors"), fig, spec, row=2, title=True)
    axs[0].set_ylabel("Compound", fontsize=16)
    add_number(axs[0], "B")

    cmap = colorby(embds, "Source")
    add_legend(fig.add_subplot(spec[4, 7]), cmap, "Source")
    axs = scatter_panel(embds.dropna(subset="colors"), fig, spec, row=3)
    axs[0].set_ylabel("Source", fontsize=16)
    add_number(axs[0], "C")

    plt.savefig(fig_path, bbox_inches="tight", dpi=300)