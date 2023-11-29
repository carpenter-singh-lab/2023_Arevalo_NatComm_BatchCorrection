include: "scib.smk"
include: "map.smk"

rule all_metrics:
    input:
        "outputs/{scenario}/metrics/{criteria}/{pipeline}_scib.parquet",
        "outputs/{scenario}/metrics/{criteria}/{pipeline}_map_negcon.parquet",
        "outputs/{scenario}/metrics/{criteria}/{pipeline}_map_nonrep.parquet",
    output:
        "outputs/{scenario}/metrics/{criteria}/{pipeline}_all_metrics.parquet",
    run:
        metrics.concat(*input, *output)

