'''Plot all figures'''
import itertools
from difflib import SequenceMatcher
import warnings

import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from plottable import ColumnDefinition, Table
from plottable.cmap import normed_cmap
from plottable.plots import bar
import seaborn as sns
from sklearn.preprocessing import minmax_scale

from .colors import BATCH_CMAP, SOURCE_CMAP, METHOD_FMT, METRIC_FMT, POSCONS, COMPOUND_COLORS, POSCON_MAP
from .ranker import Ranker


def _common_prefix_suffix(strings: list[str]):
    prefix = strings[0]
    suffix = prefix[::-1]
    for string in strings[1:]:
        match = SequenceMatcher(a=prefix, b=string).get_matching_blocks()[0]
        prefix = prefix[:match.size]
        match = SequenceMatcher(a=suffix,
                                b=string[::-1]).get_matching_blocks()[0]
        suffix = suffix[:match.size]
    suffix = suffix[::-1]
    return prefix, suffix


def _jitter(trace, win_size=0.02):
    range_x = max(trace.x) - min(trace.x)
    factor = range_x * win_size
    trace.x += np.random.uniform(-factor, factor, [len(trace.x)])

    range_y = max(trace.y) - min(trace.y)
    factor = range_y * win_size
    trace.y += np.random.uniform(-factor, factor, [len(trace.y)])


def load_all_parquet(files, key_name='file_id', placeholder='baseline'):
    dframe = []
    dframes = [pd.read_parquet(file) for file in files]
    prefix, suffix = _common_prefix_suffix(files)
    start = len(prefix)
    end = -len(suffix)
    for dframe, file in zip(dframes, files):
        dframe[key_name] = file[start:end]
    dframe = pd.concat(dframes)
    dframe[key_name] = dframe[key_name].str.strip('_').replace('', placeholder)
    return dframe


def prepare_embeddings(embd_files: list[str], output_path: str, anon=True):
    embds = load_all_parquet(embd_files, key_name='method')
    embds['method'] = embds['method'].map(lambda x: METHOD_FMT.get(x, x))

    # jitter embds
    embds_jitter = []
    for _, embd in embds.groupby('method'):
        _jitter(embd)
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=FutureWarning)
            # hack to make plotly happy in ranges
            embd[['x', 'y']] = minmax_scale(embd[['x', 'y']])
        embds_jitter.append(embd)
    embds = pd.concat(embds_jitter)

    if anon:
        # Make batches and sources annonymous
        renamer = {
            batch: f'{i:02d}'
            for i, batch in enumerate(embds['Metadata_Batch'].unique(), 1)
        }
        embds['Metadata_Batch'] = embds['Metadata_Batch'].apply(renamer.get)

        renamer = {
            source: f'{int(source.split("_")[-1]):02d}'
            for source in embds['Metadata_Source'].unique()
        }
        embds['Metadata_Source'] = embds['Metadata_Source'].apply(renamer.get)

    embds = embds.rename(
        columns={
            'Metadata_Batch': 'Batch',
            'Metadata_Source': 'Source',
            'Metadata_JCP2022': 'Compound',
            'Metadata_Row': 'Row',
            'Metadata_Column': 'Column'
        })
    embds.to_parquet(output_path, index=False)


def tidy_scores(metrics_files, metrics_redlist, methods_redlist, tidy_path):
    scores = load_all_parquet(metrics_files, key_name='method')
    scores = scores.query('metric not in @metrics_redlist')
    scores = scores[scores['method'].apply(
        lambda x: all(m not in x for m in methods_redlist))]
    scores.to_parquet(tidy_path, index=False)


def pivot_scores(tidy_path, pivot_path, micro_mean=False, macro_mean=False):
    scores = pd.read_parquet(tidy_path)
    scores['method'] = scores['method'].map(lambda x: METHOD_FMT.get(x, x))
    scores['metric'] = scores['metric'].map(lambda x: METRIC_FMT.get(x, x))
    scores = scores.pivot_table(index='method',
                                columns=['dimension', 'metric'],
                                values='score')
    scores['mean', 'batch'] = scores['batch'].mean(axis=1)
    scores['mean', 'bio'] = scores['bio'].mean(axis=1)
    if micro_mean:
        scores['mean', 'micro_mean'] = scores.mean(axis=1)
    if macro_mean:
        macro = (scores['bio'].mean(axis=1) + scores['batch'].mean(axis=1)) / 2
        scores['mean', 'macro_mean'] = macro

    # This is the default weighting from the scIB manuscript
    total = .4 * scores['mean', 'batch'] + .6 * scores['mean', 'bio']
    scores['mean', 'total'] = total
    scores = scores.sort_values(('mean', 'total'), ascending=False)
    agg_names = {
        'batch': 'Batch correction',
        'bio': 'Bio metrics',
        'total': 'Total'
    }
    scores.rename(columns=agg_names, inplace=True)
    scores.to_parquet(pivot_path)


def write_barplot(scores, title, fig_path, hue='dimension'):
    sns.reset_defaults()
    order = scores.groupby('method')['score'].mean().sort_values()[::-1]
    ax = sns.barplot(
        scores,
        x='method',
        y='score',
        hue=hue,
        order=order.index,
    )
    ax.set(title=title)
    plt.xticks(rotation=45, ha='right')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close()


def write_hbarplot(scores, title, fig_path):
    plt.figure(figsize=(6, 12))
    ax = sns.barplot(scores['mean'].reset_index().round(3),
                     y='method',
                     x='Total')
    ax.set(title=title)
    ax.bar_label(ax.containers[0], fontsize=10)
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close()


def write_umap(embds: pd.DataFrame,
               fig_path: str,
               hue,
               col_order,
               hue_order,
               palette,
               with_dmso=False):
    plt.rcParams.update({'font.size': 22})
    if not with_dmso:
        embds = embds.query('~Compound.str.startswith("DMSO")')
    g = sns.relplot(
        data=embds.sample(frac=1),  # Shuffle
        x='x',
        y='y',
        hue=hue,
        hue_order=hue_order,
        kind='scatter',
        col='method',
        col_order=col_order,
        s=6,
        col_wrap=4,
        palette=palette,
    )

    #hide axis
    g.set(xticks=[], yticks=[])
    g.set(xlabel='', ylabel='')
    sns.despine(left=True, right=True, top=True, bottom=True)

    handles, labels = g.axes[0].get_legend_handles_labels()
    g._legend.remove()
    g.figure.legend(
        handles,
        labels,
        title=f'{hue} ID',
        bbox_to_anchor=[0.95, 0.45],
        frameon=False,
        markerscale=3.5,
        ncol=1 + int(len(hue_order) > 10),
    )

    g.set_titles('{col_name}')
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close()


def barplot_all_metrics(tidy_path, fig_path):
    scores = pd.read_parquet(tidy_path)
    write_barplot(scores,
                  title='Preprocessing performance all metrics',
                  fig_path=fig_path)


def barplot_map_scores(tidy_path, fig_path):
    scores = pd.read_parquet(tidy_path)
    write_barplot(scores.query('metric.str.contains("map")'),
                  title='Preprocessing performance mAP scores',
                  fig_path=fig_path,
                  hue='metric')


def hbarplot_all_metrics(pivot_path, fig_path):
    scores = pd.read_parquet(pivot_path)
    write_hbarplot(scores, 'mean of all metrics', fig_path)


def umap_batch(embds_path, pivot_path, fig_path):
    col_order = pd.read_parquet(pivot_path).index
    embds = pd.read_parquet(embds_path)
    hue_order = np.unique(embds['Batch'])
    write_umap(embds,
               fig_path,
               'Batch',
               col_order,
               hue_order,
               BATCH_CMAP,
               with_dmso=False)


def umap_source(embds_path, pivot_path, fig_path):
    col_order = pd.read_parquet(pivot_path).index
    embds = pd.read_parquet(embds_path)
    hue_order = np.unique(embds['Source'])
    write_umap(embds,
               fig_path,
               'Source',
               col_order,
               hue_order,
               SOURCE_CMAP,
               with_dmso=False)


def umap_compound(embds_path, pivot_path, fig_path):
    col_order = pd.read_parquet(pivot_path).index
    embds = pd.read_parquet(embds_path)
    multiple_pos = embds.groupby(
        'Compound')['Metadata_Well'].nunique()[lambda x: x > 1].index
    multiple_pos = embds[embds['Compound'].isin(multiple_pos)]
    compounds = multiple_pos['Compound'].drop_duplicates().tolist()
    non_poscons = [c for c in compounds if c not in POSCONS]
    poscons = [c for c in compounds if c in POSCONS]
    hue_order = poscons + non_poscons

    compound_cmap = dict(zip(non_poscons, itertools.cycle(COMPOUND_COLORS)))
    compound_cmap.update(POSCON_MAP)
    write_umap(multiple_pos,
               fig_path,
               'Compound',
               col_order,
               hue_order,
               compound_cmap,
               with_dmso=False)


def cartesian_plane(tidy_path, min_cvar, fig_path):
    '''
    Plot scores in a cartesian plot ignoring metrics with a coefficient of variation lower than min_cvar'''
    scores = pd.read_parquet(tidy_path)
    # Find metrics with low coefficient of variation
    stats = scores.groupby('metric')['score'].agg(['std', 'mean'])
    stats['cvar'] = stats['std'] / stats['mean']
    stats = stats.fillna(0).sort_values('cvar')
    lowcvar_metrics = stats.query('cvar<@min_cvar').index.tolist()
    scores = scores.query('metric not in @lowcvar_metrics')
    Ranker.plot(scores).savefig(fig_path, bbox_inches='tight')


def best_sphering_eigen_curve(map_files, fig_path):
    comparison = []
    for map_file in map_files:
        df = pd.read_parquet(map_file)
        df.dropna(inplace=True)
        row = df['mean_average_precision'].describe()
        row['name'] = map_file[len(prefix):-len(suffix)]
        row['fraction_below_p'] = df['below_p'].sum() / len(df)
        row['fraction_below_corrected_p'] = df['below_corrected_p'].sum(
        ) / len(df)
        comparison.append(row)
    comparison = pd.DataFrame(comparison).set_index('name')
    comparison = comparison.sort_values('mean', ascending=False)
    comparison['reg'] = pd.Series(comparison.index).apply(
        lambda x: x.split('_')[-1][4:]).astype(float).values
    comparison['pipeline'] = pd.Series(
        comparison.index).apply(lambda x: x[:x.index('_reg')]).values

    ax = sns.lineplot(comparison,
                      y='mean',
                      x='reg',
                      hue='pipeline',
                      hue_order=comparison.pipeline.drop_duplicates())
    ax.set(title='Sphering lambda exploration', ylabel='mean mAP')
    plt.xscale('log')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close()


def cmap_table_fn(col_data):
    return normed_cmap(col_data, cmap=matplotlib.cm.PRGn, num_stds=2.5)


def results_table(pivot_path: str, fig_path: str):
    '''
    Adapted from:
    https://github.com/yoseflab/scib-metrics/blob/0.4.1/src/scib_metrics/benchmark/_core.py#L276-L364
    '''
    df = pd.read_parquet(pivot_path)
    column_definitions = [
        ColumnDefinition("method",
                         width=1.5,
                         textprops={
                             "ha": "left",
                             "weight": "bold"
                         }),
    ]
    score_cols = df[['Batch correction',
                     'Bio metrics']].columns.get_level_values(1)
    textprops = {"ha": "center", "bbox": {"boxstyle": "circle", "pad": 0.25}}
    groupmap = dict(df.columns.swaplevel())
    for i, col in enumerate(score_cols):
        col_def = ColumnDefinition(
            col,
            title=col.replace(" ", "\n", 1),
            textprops=textprops,
            width=1,
            cmap=cmap_table_fn(df.droplevel(0, axis='columns')[col]),
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
    agg_cols = df['mean'].columns
    for i, col in enumerate(agg_cols):
        col_def = ColumnDefinition(
            col,
            width=1,
            plot_kw=plot_kw,
            title=col.replace(" ", "\n", 1),
            plot_fn=bar,
            group='Aggregate score',
            border="left" if i == 0 else None,
        )
        column_definitions.append(col_def)

    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(len(df.columns) * 1.25, 3 + 0.3 * len(df)))
    tab = Table(
        df.droplevel(0, axis='columns').reset_index(),
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
    fig.savefig(fig_path,
                facecolor=ax.get_facecolor(),
                dpi=300,
                bbox_inches='tight')
