import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from .data import _common_prefix_suffix
from .ranker import Ranker


def cartesian_plane(tidy_path, min_cvar, fig_path):
    """
    Plot scores in a cartesian plot ignoring metrics with a coefficient of variation lower than min_cvar"""
    scores = pd.read_parquet(tidy_path)
    # Find metrics with low coefficient of variation
    stats = scores.groupby("metric")["score"].agg(["std", "mean"])
    stats["cvar"] = stats["std"] / stats["mean"]
    stats = stats.fillna(0).sort_values("cvar")
    lowcvar_metrics = stats.query("cvar<@min_cvar").index.tolist()
    scores = scores.query("metric not in @lowcvar_metrics")
    Ranker.plot(scores).savefig(fig_path, bbox_inches="tight")


def best_sphering_eigen_curve(map_files, fig_path):
    comparison = []
    (
        prefix,
        suffix,
    ) = _common_prefix_suffix(map_files)
    for map_file in map_files:
        df = pd.read_parquet(map_file)
        df.dropna(inplace=True)
        row = df["mean_average_precision"].describe()
        row["name"] = map_file[len(prefix) : -len(suffix)]
        row["fraction_below_p"] = df["below_p"].sum() / len(df)
        row["fraction_below_corrected_p"] = df["below_corrected_p"].sum() / len(df)
        comparison.append(row)
    comparison = pd.DataFrame(comparison).set_index("name")
    comparison = comparison.sort_values("mean", ascending=False)
    comparison["reg"] = (
        pd.Series(comparison.index)
        .apply(lambda x: x.split("_")[-1][4:])
        .astype(float)
        .values
    )
    comparison["pipeline"] = (
        pd.Series(comparison.index).apply(lambda x: x[: x.index("_reg")]).values
    )

    ax = sns.lineplot(
        comparison,
        y="mean",
        x="reg",
        hue="pipeline",
        hue_order=comparison.pipeline.drop_duplicates(),
    )
    ax.set(title="Sphering lambda exploration", ylabel="mean mAP")
    plt.xscale("log")
    plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    plt.savefig(fig_path, bbox_inches="tight")
    plt.close()
