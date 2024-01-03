'''Plot all figures'''
from difflib import SequenceMatcher
import warnings

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import minmax_scale

from .colors import BATCH_CMAP, SOURCE_CMAP
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


def pivot_scores(tidy_path, pivot_path):
    scores = pd.read_parquet(tidy_path)
    scores = scores.pivot_table(index='method',
                                columns=['dimension', 'metric'],
                                values='score')
    scores['mean', 'batch'] = scores['batch'].mean(axis=1)
    scores['mean', 'bio'] = scores['bio'].mean(axis=1)
    scores['mean', 'micro_mean'] = scores.mean(axis=1)
    scores['mean', 'macro_mean'] = (scores['bio'].mean(axis=1) +
                                    scores['batch'].mean(axis=1)) / 2
    scores = scores.sort_values(('mean', 'macro_mean'), ascending=False)
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
                     x='macro_mean')
    ax.set(title=title)
    ax.bar_label(ax.containers[0], fontsize=10)
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close()


def write_umap(embds_path, fig_path, hue, order, palette, with_dmso=False):
    embds = pd.read_parquet(embds_path)
    plt.rcParams.update({'font.size': 22})
    if not with_dmso:
        embds = embds.query('~Compound.str.startswith("DMSO")')
    g = sns.relplot(
        data=embds.sample(frac=1),  # Shuffle
        x='x',
        y='y',
        hue=hue,
        hue_order=embds[hue].drop_duplicates().sort_values(),
        kind='scatter',
        col='method',
        col_order=order,
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
    write_umap(embds_path,
               fig_path,
               'Batch',
               col_order,
               BATCH_CMAP,
               with_dmso=False)


def umap_source(embds_path, pivot_path, fig_path):
    col_order = pd.read_parquet(pivot_path).index
    write_umap(embds_path,
               fig_path,
               'Source',
               col_order,
               SOURCE_CMAP,
               with_dmso=False)


def cartesian_plane(tidy_path, fig_path):
    scores = pd.read_parquet(tidy_path)
    Ranker.plot(scores).savefig(fig_path)


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
