import pandas as pd

from metrics import DIMENSION_MAP
from plot.colors import METRIC_FMT


def get_table():
    path = "./outputs/scenario_{}/metrics/{}/mad_int_featselect_bbknn_scib.parquet"
    scores = []
    for scn in range(1, 6):
        if scn in (3, 5):
            criteria = "prod"
        else:
            criteria = "target2"
        scib_path = path.format(scn, criteria)
        df = pd.read_parquet(scib_path).set_index("metric")["score"]
        df = df.reset_index()
        df["dimension"] = df["metric"].map(DIMENSION_MAP)
        df["scenario"] = f"Scenario {scn}"
        scores.append(df)
    df = pd.concat(scores).reset_index(drop=True)
    df["metric"] = df["metric"].map(lambda x: METRIC_FMT.get(x, x))
    df = df.pivot_table(
        index="scenario", columns=["dimension", "metric"], values="score"
    ).round(3)
    df.columns = df.columns.set_levels(["Batch correction", "Bio metrics"], level=0)
    return df


df = get_table()
