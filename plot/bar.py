import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def write_barplot(scores, title, fig_path, hue="dimension"):
    sns.reset_defaults()
    order = scores.groupby("method")["score"].mean().sort_values()[::-1]
    ax = sns.barplot(
        scores,
        x="method",
        y="score",
        hue=hue,
        order=order.index,
    )
    ax.set(title=title)
    plt.xticks(rotation=45, ha="right")
    plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    plt.savefig(fig_path, bbox_inches="tight")
    plt.close()


def write_hbarplot(scores, title, fig_path):
    plt.figure(figsize=(6, 12))
    ax = sns.barplot(scores["mean"].reset_index().round(3),
                     y="method",
                     x="Overall")
    ax.set(title=title)
    ax.bar_label(ax.containers[0], fontsize=10)
    plt.savefig(fig_path, bbox_inches="tight")
    plt.close()


def all_metrics(tidy_path, fig_path):
    scores = pd.read_parquet(tidy_path)
    write_barplot(scores,
                  title="Preprocessing performance all metrics",
                  fig_path=fig_path)


def map_scores(tidy_path, fig_path):
    scores = pd.read_parquet(tidy_path)
    write_barplot(
        scores.query('metric.str.contains("map")'),
        title="Preprocessing performance mAP scores",
        fig_path=fig_path,
        hue="metric",
    )


def all_metrics_h(pivot_path, fig_path):
    scores = pd.read_parquet(pivot_path)
    write_hbarplot(scores, "mean of all metrics", fig_path)
