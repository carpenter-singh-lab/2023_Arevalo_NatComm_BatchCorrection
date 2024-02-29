import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

methods = [
        "mad_drop_int_featselect",
        "mad_int_featselect",
        "mad_drop_int",
        "mad_drop_featselect",
        "mad_int",
        "mad_drop",
        "mad_featselect",
        "mad_clip",
        "mad_imputeknn",
        "mad_imputemedian",
        "mad"
        ]

df = []
for method in methods:
    sdf = pd.read_parquet(f'./outputs/scenario_1/metrics/target2/{method}_all_metrics.parquet')
    sdf['method'] = method
    df.append(sdf)
df = pd.concat(df)

df = df.query('metric.str.contains("map")')
ax = sns.barplot(df, x="method", y="score", hue="metric")
ax.set_xlabel(None)
ax.set(title="Preprocessing mAP scores")
plt.xticks(rotation=45, ha="right")
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.savefig("figures/sup_figure_A.pdf", bbox_inches='tight')
