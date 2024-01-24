import base64
from io import BytesIO

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def show_inline():
    iobytes = BytesIO()
    plt.savefig(iobytes, format="png", bbox_inches="tight", dpi=300)
    iobytes.seek(0)
    encode = base64.b64encode(iobytes.read())
    props = "File=inline=1"
    print(f"\x1b]1337;{props}:" + encode.decode() + "\x07\n")


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


def add_table(ax: plt.Axes):
    columns = [f"Scenario {i}" for i in scenarios]
    rows = ["Sources", "Microscopes", "Compounds"]
    cell_text = [[1, 3, 3, 5, 5], [1, 1, 1, 3, 3],
                 [306, 306, "80000+", 306, "80000+"]]
    ax.table(cellText=cell_text,
             rowLabels=rows,
             colLabels=columns,
             loc="bottom")

    ax.set(xlabel=None)
    ax.xaxis.set_ticks([])
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
    ax.set_title("Performance distribution - All methods")
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
