import itertools

import pandas as pd
from matplotlib import pyplot as plt


def despine(ax: plt.Axes):
    for side in ["left", "right", "top", "bottom"]:
        ax.spines[side].set_visible(False)
        ax.spines[side].set_visible(False)
    ax.tick_params(which="both",
                   bottom=False,
                   left=False,
                   labelbottom=False,
                   labelleft=False)


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
