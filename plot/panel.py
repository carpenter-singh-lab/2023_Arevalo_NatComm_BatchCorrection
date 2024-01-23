from functools import partial

import pandas as pd
from matplotlib import pyplot as plt

from .colors import MICRO_CMAP, SOURCE_CMAP
from .data import query_multiple_pos
from .scatter import create_compound_cmap, despine, scatter_panel
from .table import add_colorbar
from .table import draw as draw_table


def full(embd_path: str, pivot_path: str, fig_path: str):
    embds = pd.read_parquet(embd_path)
    embds = embds.query('~Compound.str.startswith("DMSO")')

    # Sort based on the pivot scores table
    rank = pd.read_parquet(pivot_path).index.to_list()
    embds = embds.set_index("Method").loc[rank].reset_index()

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
    cpd_order = [c.replace('JCP2022_', '') for c in cpd_order]
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
    plt.savefig(fig_path, bbox_inches="tight")
