import itertools

import pandas as pd
from matplotlib import pyplot as plt

from plot.colors import COMPOUND_COLORS, POSCON_MAP, POSCONS


def despine(ax: plt.Axes):
    for side in ["left", "right", "top", "bottom"]:
        ax.spines[side].set_visible(False)
        ax.spines[side].set_visible(False)
    ax.tick_params(which="both",
                   bottom=False,
                   left=False,
                   labelbottom=False,
                   labelleft=False)


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
    methods = embds["Method"].drop_duplicates().to_list()
    for i, method in enumerate(methods):
        points = embds.query("Method==@method").sample(frac=1)
        x, y = points["x"], points["y"]
        colors = points["colors"]
        ax = fig.add_subplot(spec[row, i])
        ax.scatter(x, y, c=colors, s=6)
        despine(ax)
        if title:
            ax.set_title(method)
