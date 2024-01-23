import matplotlib
import pandas as pd
from matplotlib import pyplot as plt
from plottable import ColumnDefinition, Table
from plottable.plots import bar


def add_colorbar(fig, ax_colorbar):
    fig.colorbar(
        get_scalar_mapppable([]),
        cax=ax_colorbar,
        orientation="horizontal",
    )
    ax_colorbar.xaxis.set_ticks_position("top")
    ax_colorbar.set_ylabel("Score ", ha="right", rotation="horizontal")


def white_yellow_green_cm():
    lut_size = 256
    spec = [
        (1.0, 1.0, 1.0),
        (0.90196078431372551, 0.96078431372549022, 0.81568627450980391),
        (0.30196078431372547, 0.5725490196078431, 0.12941176470588237),
    ]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("wYlGn", spec, lut_size)
    return cmap


def get_scalar_mapppable(col_data, norm_type=None):
    if norm_type == "minmax":
        vmin = col_data.min()
        vmax = col_data.max()
    if norm_type == "interquartile":
        # taken from plottable.cmap.normed_cmap
        num_stds = 2.5
        _median, _std = col_data.median(), col_data.std()
        vmin = _median - num_stds * _std
        vmax = _median + num_stds * _std
    else:
        vmin, vmax = 0, 1

    cmap = white_yellow_green_cm()
    norm = matplotlib.colors.Normalize(vmin, vmax)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    return m


def draw(pivot_path: str, ax: plt.Axes):
    """
    Adapted from:
    https://github.com/yoseflab/scib-metrics/blob/0.4.1/src/scib_metrics/benchmark/_core.py#L276-L364
    """
    df = pd.read_parquet(pivot_path)

    column_definitions = [
        ColumnDefinition(
            "method", width=1.5, textprops={"ha": "left", "weight": "bold"}
        ),
    ]
    score_cols = df[["Batch correction", "Bio metrics"]].columns.get_level_values(1)
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
        textprops={"fontsize": 10, "ha": "center"},
        row_divider_kw={"linewidth": 1, "linestyle": (0, (1, 5))},
        col_label_divider_kw={"linewidth": 1, "linestyle": "-"},
        column_border_kw={"linewidth": 1, "linestyle": "-"},
        index_col="method",
    )
    tab.autoset_fontcolors(colnames=list(df.columns.get_level_values(1)))
