from plot import proj


rule mde:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/projection/{pipeline}_mde.parquet",
    run:
        proj.mde(*input, *output)


rule pca:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/projection/{pipeline}_pca.parquet",
    run:
        proj.pca(*input, *output)


rule umap:
    input:
        f"outputs/{{scenario}}/metrics/{criteria}/scib/{{pipeline}}_clusters.h5ad",
    output:
        "outputs/{scenario}/projection/{pipeline}_umap.parquet",
    run:
        proj.umap(*input, *output)
