from difflib import SequenceMatcher
from glob import glob

from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns

import plot


def plot_best_sphering(map_files, fig_path):
    comparison = []
    for map_file in map_files:
        df = pd.read_parquet(map_file)
        df.dropna(inplace=True)
        row = df['mean_average_precision'].describe()
        row['name'] = map_file[len(prefix):-len(suffix)]
        row['fraction_below_p'] = df['below_p'].sum() / len(df)
        row['fraction_below_corrected_p'] = df['below_corrected_p'].sum() / len(df)
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
    ax.set(title=f'Sphering lambda exploration', ylabel='mean mAP')
    plt.xscale('log')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close()


def common_prefix_suffix(strings: list[str]):
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


def load_all_parquet(files, key_name='file_id'):
    dframe = []
    dframes = [pd.read_parquet(file) for file in files]
    prefix, suffix = common_prefix_suffix(files)
    start = len(prefix)
    end = -len(suffix)
    for dframe, file in zip(dframes, files):
        dframe[key_name] = file[start:end]
    dframe = pd.concat(dframes)
    return dframe


def load_scores(metrics_files, metrics_redlist, methods_redlist):
    scores = load_all_parquet(metrics_files, key_name='method')
    scores = scores.query('metric not in @metrics_redlist')
    scores = scores[scores['method'].apply(
        lambda x: all(m not in x for m in methods_redlist))]
    return scores


def pivot(scores):
    scores = scores.pivot_table(index='method',
                                columns=['dimension', 'metric'],
                                values='score')
    scores['mean', 'batch'] = scores['batch'].mean(axis=1)
    scores['mean', 'bio'] = scores['bio'].mean(axis=1)
    scores['mean', 'micro_mean'] = scores.mean(axis=1)
    scores['mean', 'macro_mean'] = (scores['bio'].mean(axis=1) + scores['batch'].mean(axis=1)) / 2
    scores = scores.sort_values(('mean', 'macro_mean'), ascending=False)
    return scores


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
    ax = sns.barplot(scores['mean'].reset_index().round(3), y='method', x='macro_mean')
    ax.set(title=title)
    ax.bar_label(ax.containers[0], fontsize=10)
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close()


def write_umap(embd_files, fig_path, hue, order, palette, with_dmso=False):
    embds = plot.embeddings.load_embeddings(embd_files)
    prefix, suffix = common_prefix_suffix(embd_files)
    embds['method'] = embds['path'].str[len(prefix):-len(suffix)]
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



redlist = [
    'pcr', 'pcr_batch', 'il_f1', 'il_asw', 'negcon_fraction_below_p',
    'negcon_fraction_below_corrected_p', 'nonrep_fraction_below_p',
    'nonrep_fraction_below_corrected_p', 'lisi_label'
],
corr_meth = ['harmony', 'scanorama', 'mnn', 'combat']
corr_meth = []
scenario = 'scenario_1'

metrics_files = glob(f'outputs/{scenario}/metrics/target2/*_all_metrics.parquet')
prefix, suffix = common_prefix_suffix(metrics_files)
scores = load_scores(metrics_files, redlist, corr_meth)
write_barplot(
    scores,
    title='Preprocessing performance all metrics',
    fig_path=f'{scenario}_all_metrics_barplot.png',
)

write_barplot(
    scores.query('metric.str.contains("map")'),
    title='Preprocessing performance mAP scores',
    fig_path=f'{scenario}_map_scores_barplot.png',
    hue='metric'
)

pivot_scores = pivot(scores)
write_hbarplot(pivot_scores, 'mean of all metrics', f'{scenario}_mean_all_metrics.png')

embd_files = glob(f'outputs/{scenario}/projection/*_umap.parquet')
col_order = pivot_scores.index.str[:-1].values
col_order = [c if len(c) else 'baseline' for c in col_order]

write_umap(embd_files, f'{scenario}_source_umap.png', 'Source', col_order, plot.SOURCE_CMAP, with_dmso=False)
write_umap(embd_files, f'{scenario}_batch_umap.png', 'Batch', col_order, plot.BATCH_CMAP, with_dmso=False)
write_umap(embd_files, f'{scenario}_compound_umap.png', 'Compound', col_order, None, with_dmso=False)
