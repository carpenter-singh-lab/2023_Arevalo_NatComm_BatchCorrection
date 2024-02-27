import base64
from io import BytesIO

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def show_inline(close=False):
    iobytes = BytesIO()
    plt.savefig(iobytes, format="png", bbox_inches="tight", dpi=300)
    iobytes.seek(0)
    encode = base64.b64encode(iobytes.read())
    props = "File=inline=1"
    print(f"\x1b]1337;{props}:" + encode.decode() + "\x07\n")
    if close:
        plt.close()


scenarios = [1, 2, 3, 4, 5]
nsources = dict(zip(scenarios, [1, 3, 3, 5, 5]))
nmicros =dict(zip(scenarios, [1, 1, 1, 3, 3]))
ncomps = dict(zip(scenarios, [306, 306, "80000+", 306, "80000+"]))
# swap order to make it increasing difficulty
scenarios = [1, 2, 4, 3, 5]

scores = []
min_cvar = 0.1
lowcvar_metrics = []
for scn in [f"scenario_{i}" for i in scenarios]:
    df = pd.read_parquet(f"outputs/{scn}/plots/data/tidy_scores.parquet")
    # Find metrics with low coefficient of variation
    stats = df.groupby("metric")["score"].agg(["std", "mean"])
    stats["cvar"] = stats["std"] / stats["mean"]
    stats = stats.fillna(0).sort_values("cvar")
    lowcvar_metrics += stats.query("cvar<@min_cvar").index.tolist()
    df["scenario"] = scn
    scores.append(df)
scores = pd.concat(scores)
# scores = scores.query("metric not in @lowcvar_metrics")
scores["rank"] = scores.groupby(["scenario",
                                 "metric"])["score"].rank().astype(int)
scores["dimension"] = scores["dimension"].apply(lambda x: {
    "batch": "Batch correction",
    "bio": "Bio metrics"
}.get(x, x))

box_scores = scores  # .query("method=='harmony'")


def add_table(ax: plt.Axes, add_legend=True, row_labels=True):
    columns = [f"Scenario {i}" for i in scenarios]
    rows = ["Sources", "Microscopes", "Compounds"]
    cell_text = [[nsources[i] for i in scenarios],
                 [nmicros[i] for i in scenarios],
                 [ncomps[i] for i in scenarios],
                ]

    tbl = ax.table(cellText=cell_text,
             rowLabels=rows if row_labels else None,
             colLabels=columns,
             loc="bottom")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)

    ax.set(xlabel=None)
    ax.xaxis.set_ticks([])
    if add_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            handles[:2],
            labels[:2],
            title="Dimension",
            bbox_to_anchor=(1, 1),
            loc="upper left",
        )


def plot_box_scores(box_scores: pd.DataFrame):
    ax = sns.boxplot(
        box_scores,
        hue="dimension",
        y="score",
        x="scenario",
        gap=0.3,
        fill=False,
    )
    ax.set_title("Harmony improvement over the Baseline")
    sns.stripplot(
        box_scores,
        hue="dimension",
        y="score",
        alpha=0.4,
        x="scenario",
        dodge=True,
        ax=ax,
    )
    add_table(ax)
    ax.set_ylabel('Harmony_score - Baseline_score')
    show_inline()
    plt.close()


def plot_metrics_agg():
    facet_col = "dimension"
    facet_cols = scores[facet_col].unique()
    fig, axes = plt.subplots(nrows=1, ncols=len(facet_cols))
    scores["scenario"] = scores["scenario"].apply(lambda x: {
        "scenario_3": "scenario_4",
        "scenario_4": "scenario_3"
    }.get(x, x))
    for i, col in enumerate(facet_cols):
        sns.barplot(
            scores.query(f"{facet_col}==@col").sort_values("scenario"),
            x="metric",
            y="score",
            hue="scenario",
            legend=i == 1,
            ax=axes[i],
        ).set_title(col)
    axes[-1].legend(loc="upper right", bbox_to_anchor=(-0.05, 1))
    show_inline()
    plt.close()

    facet_col = "dimension"
    facet_cols = scores[facet_col].unique()
    fig, axes = plt.subplots(nrows=1, ncols=len(facet_cols))
    scores["scenario"] = scores["scenario"].apply(lambda x: {
        "scenario_3": "scenario_4",
        "scenario_4": "scenario_3"
    }.get(x, x))
    for i, col in enumerate(facet_cols):
        sns.barplot(
            scores.query(f"{facet_col}==@col").sort_values("scenario"),
            x="scenario",
            y="rank",
            hue="method",
            legend=i == 1,
            ax=axes[i],
        )
    axes[-1].legend(loc="upper right", bbox_to_anchor=(-0.05, 1))
    show_inline()
    plt.close()

fig = plt.figure(figsize=(10, 4))
spec = fig.add_gridspec(1, 2)
delta = (scores.query('method=="harmony"').set_index(['scenario', 'metric',
                                                   'dimension']).score -
         scores.query('method=="baseline"').set_index(['scenario', 'metric',
                                                       'dimension']).score).reset_index()
delta2 = (scores.query('method=="desc"').set_index(['scenario', 'metric',
                                                       'dimension']).score -
          scores.query('method=="baseline"').set_index(['scenario', 'metric',
                                                        'dimension']).score).reset_index()

ax = fig.add_subplot(spec[0, 0])
sns.boxplot(delta, x='scenario', y='score', hue='dimension', fill=False, ax=ax)
sns.stripplot(delta, x='scenario', y='score', hue='dimension', dodge=True, ax=ax, alpha=0.4)
add_table(ax, add_legend=False)
ax.set_ylabel("Improvement over Baseline")
# ax.get_legend().remove()
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], labels[:2], ncol=2)
ax.set_title('Harmony (best)')
ax = fig.add_subplot(spec[0, 1], sharey=ax)
plt.setp(ax.get_yticklabels(), visible=False)
sns.boxplot(delta2, x='scenario', y='score', hue='dimension', fill=False, ax=ax)
sns.stripplot(delta2, x='scenario', y='score', hue='dimension', dodge=True, ax=ax, alpha=0.4)
add_table(ax, add_legend=False, row_labels=False)
ax.set_title('DESC (worst)')
ax.get_legend().remove()
fig.tight_layout()
show_inline()
plt.close('all')
for scn, data in scores.groupby('scenario'):
    agg_scores = data.groupby('dimension')['score'].agg(['mean', 'std'])
    plt.errorbar(x=agg_scores.iloc[0, 0], y=agg_scores.iloc[1, 0], xerr=agg_scores.iloc[0, 1], yerr=agg_scores.iloc[1, 1], label=scn)
plt.legend()
plt.gca().set_ylabel('Bio metrics')
plt.gca().set_xlabel('Batch correction')
show_inline()
plt.close('all')
