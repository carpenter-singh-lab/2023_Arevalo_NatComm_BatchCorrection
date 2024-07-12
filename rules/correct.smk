include: "sphering.smk"


rule combat:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_combat.parquet",
    params:
        batch_key=lambda config: config["batch_key"] if isinstance(config["batch_key"], string) else config["batch_key"][0],
    run:
        correct.combat(*input, *params, *output)


rule harmony:
    input:
        data="outputs/{scenario}/{pipeline}.parquet",
        script="scripts/harmony.py"
    output:
        "outputs/{scenario}/{pipeline}_harmony.parquet"
    log:
        "logs/{scenario}/{pipeline}_harmony.log"
    conda:
        "../envs/harmony.yaml"  
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0]
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && python {input.script} harmony {input.data} {params.batch_key} {output} &> {log}
        """


rule pca_harmony:
    input:
        data="outputs/{scenario}/{pipeline}.parquet",
        script="scripts/harmony.py"
    output:
        "outputs/{scenario}/{pipeline}_pca_harmony.parquet"
    log:
        "logs/{scenario}/{pipeline}_pca_harmony.log"
    conda:
        "../envs/harmony.yaml"  
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0]
    shell:
        "python {input.script} pca_harmony {input.data} {params.batch_key} {output}"


rule scanorama:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_scanorama.parquet",
    params:
        batch_key=lambda config: config["batch_key"] if isinstance(config["batch_key"], string) else config["batch_key"][0],
    run:
        correct.scanorama(*input, *params, *output)


rule pca_scanorama:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_pca_scanorama.parquet",
    params:
        batch_key=lambda config: config["batch_key"] if isinstance(config["batch_key"], string) else config["batch_key"][0],
    run:
        correct.pca_scanorama(*input, *params, *output)


rule mnn:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_mnn.parquet",
    params:
        batch_key=lambda config: config["batch_key"] if isinstance(config["batch_key"], string) else config["batch_key"][0],
    run:
        correct.mnn(*input, *params, *output)


rule desc:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_desc.parquet",
    params:
        batch_key=lambda config: config["batch_key"] if isinstance(config["batch_key"], string) else config["batch_key"][0],
    run:
        correct.desc(*input, *params, *output)

rule scvi:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_scvi.parquet",
    log:
        "{scenario}/{pipeline}_scvi.log"
    conda:
        "./../envs/scvi.yaml"
    params:
        batch_key=config["batch_key"],
        label_key=config["label_key"],
    script:
        "scripts/scvi.py {input} {params.batch_key} {params.label_key} {output}"

rule sysvi:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_sysvi.parquet",
    params:
        batch_key=config["batch_key"],
        label_key=config["label_key"],
    run:
        correct.sysvi(*input, *params, *output)


rule fastMNN:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_fastMNN.parquet",
    params:
        batch_key=lambda config: config["batch_key"] if isinstance(config["batch_key"], string) else config["batch_key"][0],
    shell:
        "Rscript correct/fastMNN.R {input} {output} {params.batch_key}"


rule seurat:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_seurat_{seurat_method}.parquet",
    params:
        batch_key=lambda config: config["batch_key"] if isinstance(config["batch_key"], string) else config["batch_key"][0],
    shell:
        "Rscript correct/seurat.R {input} {output} {params.batch_key} {wildcards.seurat_method}"
