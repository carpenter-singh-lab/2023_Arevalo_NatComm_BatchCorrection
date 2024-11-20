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

    # prepare cmaps so we can pass them as callables
    Blues = matplotlib.cm.get_cmap("Blues")
    Greens = matplotlib.cm.get_cmap("Greens")

    # extract labels from columns so we can construct queries later
    print(df)
    eval_keys = df.columns.get_level_values("eval_key").unique()
    metric_type_keys = df.columns.get_level_values("metric_type").unique()
    metric_keys = df.columns.get_level_values("metric").unique().drop(["mean_batch", "mean_bio", "mean_overall"])

    # flatten the MultiIndex so we can properly create a plottable
    df_for_table = df.reset_index()
    df_for_table.columns = [
        "_".join(filter(None, map(str, col))).strip() if isinstance(col, tuple) else col
        for col in df_for_table.columns.values
    ]
    print(df_for_table.columns)
    # df_for_table.fillna(0, inplace=True)  # TMP HACK

    # Add first column for the method names
    column_definitions = [
        ColumnDefinition(
            "Method", width=1.5, textprops={"ha": "left", "weight": "bold"}
        ),
    ]

    # define how circles and barplots look
    circle_props = {"ha": "center", "bbox": {"boxstyle": "circle", "pad": 0.25}}
    barplot_props = {
        "cmap": matplotlib.cm.YlGnBu,
        "plot_bg_bar": False,
        "annotate": True,
        "height": 0.9,
        "formatter": "{:.2f}",
    }

    # we'll use this to resort the df columns later
    final_col_order = []

    # iterate over keys and generate columns for each
    for eval_key in eval_keys:
        # we'll iterate over this later for the aggregated scores
        mean_col_fstring_tuples = []

        for metric_type in metric_type_keys:
            print(metric_type)
            # Skip 'aggregate_score' as it contains mean values
            if metric_type == 'aggregate_score':
                continue

            for metric in metric_keys:
                colname = f"{eval_key}_{metric_type}_{metric}"

                # we might not have all metrics for all eval_keys (esp during testing)
                if colname in df_for_table.columns:
                    col_def = ColumnDefinition(
                        colname,
                        title=metric.replace("_", "\n", 1),
                        textprops=circle_props,
                        width=1,
                        cmap=Blues if metric_type == "batch_correction" else Greens,
                        group=eval_key,
                        formatter="{:.2f}",
                    )
                    column_definitions.append(col_def)
                    final_col_order.append(colname)

            # Collect mean columns
            metric_substring = metric_type.split("_")[0]
            mean_colname = f"{eval_key}_aggregate_score_mean_{metric_substring}"
            if mean_colname in df_for_table.columns:
                mean_col_fstring_tuples.append((eval_key, 'aggregate_score', f"mean_{metric_substring}"))

        # Now, append 'mean_overall'
        mean_overall_colname = f"{eval_key}_aggregate_score_mean_overall"
        if mean_overall_colname in df_for_table.columns:
            mean_col_fstring_tuples.append((eval_key, 'aggregate_score', "mean_overall"))

        # Add mean columns for reordering
        # final_col_order.extend([f"{ek}_{mt}_{meancol}" for ek, mt, meancol in mean_col_fstring_tuples])

        title_lookup = {
            "mean_batch": "mean\nbatch",
            "mean_bio": "mean\nbio",
            "mean_overall": f"mean\n{eval_key.replace('Metadata_', '')}",
        }

        border_lookup = {
            0: "left",
            len(mean_col_fstring_tuples) - 1: "right",
        }

        for i, (ek, mt, meancol) in enumerate(mean_col_fstring_tuples):
            colname = f"{ek}_{mt}_{meancol}"
            col_def = ColumnDefinition(
                colname,
                title=title_lookup.get(meancol, meancol).replace("_", "\n", 1),
                plot_kw=barplot_props,
                width=1,
                plot_fn=bar,
                group=eval_key,
                formatter="{:.2f}",
                border=border_lookup.get(i),
            )
            column_definitions.append(col_def)
            final_col_order.append(colname)

    # Calculate mean-score across all eval_keys and add it to the table
    name_for_col_that_averages_across_eval_keys = "total_mean"
    total_mean = df_for_table[[col for col in df_for_table.columns if "aggregate_score_mean_overall" in col]].mean(axis=1)
    df_for_table["total_mean"] = total_mean

    col_def = ColumnDefinition(
        name_for_col_that_averages_across_eval_keys,
        title="Total",
        plot_kw=barplot_props,
        width=1,
        plot_fn=bar,
        formatter="{:.2f}",
        border="both",
    )
    column_definitions.append(col_def)
    final_col_order.append(name_for_col_that_averages_across_eval_keys)

    # Reorder columns and rows
    df_for_table = df_for_table[["Method"] + final_col_order]
    df_for_table = df_for_table.sort_values(by=name_for_col_that_averages_across_eval_keys, ascending=False)

    # create the table
    plt.style.use("default")
    tab = Table(
        df_for_table,
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
        index_col="Method",
    )
    # Adjust colnames for font color settings
    colnames = [col_def.name for col_def in column_definitions[1:]]  # Exclude "Method"
    tab.autoset_fontcolors(colnames=colnames)