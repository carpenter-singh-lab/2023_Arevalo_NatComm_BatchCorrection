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
        adata_path="outputs/{scenario}/{pipeline}_all_methods.h5ad",
    output:
        "outputs/{scenario}/projection/{pipeline}_umap.parquet",
        expand("outputs/{{scenario}}/projection/{{pipeline}}_{method}_umap.parquet", method=METHODS),
    run:
        proj.umap(*input, *output)
