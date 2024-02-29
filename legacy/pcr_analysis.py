from glob import glob

import numpy as np
import pandas as pd

from plot.data import _common_prefix_suffix

df = []
for scn in range(1, 6):
    criteria = "target2" if scn in (1, 2, 4) else "prod"
    paths = glob(f"outputs/scenario_{scn}/metrics/{criteria}/scib/*pcr*")
    pre, suf = _common_prefix_suffix(paths)
    names = [p[len(pre) : -len(suf)] for p in paths]
    values = [np.fromfile(c).item() for c in paths]
    scores = pd.DataFrame({"name": names, "value": values})
    scores["scenario"] = scn
    df.append(scores)
df = pd.concat(df)
df["method"] = df["name"].str.split("_", n=1, expand=True)[0]
df["metric"] = df["name"].str.split("_", n=1, expand=True)[1]
df["metric"] = df["metric"].replace(["batch"], "pcr_batch")
df["metric"] = df["metric"].fillna("pcr")
df["method"] = df["method"].replace(["pcr"], "mad", regex=True)

df.loc[df["metric"] == "pcr", "value"] = 1 - df.loc[df["metric"] == "pcr", "value"]
df.loc[(df["method"] == "mad") & (df["metric"] == "pcr_batch"), "value"] = None
pivot = df.pivot_table(index="method", columns=["scenario", "metric"], values="value")
pivot.to_csv('pcr_scores.csv')
