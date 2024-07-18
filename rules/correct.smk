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
        script="scripts/correct_with_harmony.py"
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
        script="scripts/correct_with_harmony.py"
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
        data="outputs/{scenario}/{pipeline}.parquet",
        script="scripts/correct_with_scanorama.py"
    output:
        "outputs/{scenario}/{pipeline}_scanorama.parquet"
    log:
        "logs/{scenario}/{pipeline}_scanorama.log"
    conda:
        "../envs/scanorama.yaml"  
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0]
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && python {input.script} scanorama {input.data} {params.batch_key} {output} &> {log}
        """



rule pca_scanorama:
    input:
        data="outputs/{scenario}/{pipeline}.parquet",
        script="scripts/correct_with_scanorama.py"
    output:
        "outputs/{scenario}/{pipeline}_pca_scanorama.parquet"
    log:
        "logs/{scenario}/{pipeline}_pca_scanorama.log"
    conda:
        "../envs/scanorama.yaml"  
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0]
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && python {input.script} pca_scanorama {input.data} {params.batch_key} {output} &> {log}
        """


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
        data="outputs/{scenario}/{pipeline}.parquet",
        script="scripts/correct_with_scvi.py"
    output:
        "outputs/{scenario}/{pipeline}_scvi.parquet",
    log:
        "logs/{scenario}/{pipeline}_scvi.log"
    conda:
        "../envs/scvi.yaml"  
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        label_key=config["label_key"],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && python {input.script} {input.data} {params.batch_key} {params.label_key} {output} &> {log}
        """


rule sysvi:
    input:
        data="outputs/{scenario}/{pipeline}.parquet",
        script="scripts/correct_with_sysvi.py"
    output:
        "outputs/{scenario}/{pipeline}_sysvi.parquet",
    log:
        "logs/{scenario}/{pipeline}_sysvi.log"
    conda:
        "../envs/sysvi.yaml"  
    params:
        batch_key=config["batch_key"],
        label_key=config["label_key"],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && python {input.script} {input.data} '{params.batch_key}' {params.label_key} {output} &> {log}
        """

rule scpoli:
    input:
        data="outputs/{scenario}/{pipeline}.parquet",
        script="scripts/correct_with_scpoli.py"
    output:
        "outputs/{scenario}/{pipeline}_scpoli.parquet",
    log:
        "logs/{scenario}/{pipeline}_scpoli.log"
    conda:
        "../envs/scpoli.yaml"  
    params:
        batch_key=config["batch_key"],
        label_key=config["label_key"],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && python {input.script} {input.data} '{params.batch_key}' {params.label_key} {output} &> {log}
        """


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
        data="outputs/{scenario}/{pipeline}.parquet",
        script="scripts/correct_with_seurat.R",
        method=expand("{method}", method=config["seurat_methods"])
    output:
        "outputs/{scenario}/{pipeline}_seurat_{input.method}.parquet",
    log:
        "logs/{scenario}/{pipeline}_seurat_{input.method}.log"
    conda:
        "../envs/seurat.yaml"  
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && Rscript {input.script} {input.data} {params.batch_key} {input.method} {output} &> {log}
        """
