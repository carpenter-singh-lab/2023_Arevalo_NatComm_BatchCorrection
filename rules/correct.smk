include: "sphering.smk" # not modifying this for now, dependencies are unproblematic

rule combat:
    input:
        data="outputs/{scenario}/{pipeline}.parquet",
        script="scripts/correct_with_combat.py"
    output:
        "outputs/{scenario}/{pipeline}_combat.parquet"
    log:
        "logs/{scenario}/{pipeline}_combat.log"
    conda:
        "../envs/harmony.yaml"  # we only need scanpy so this will do
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && python {input.script} {input.data} {params.batch_key} {output} &> {log}
        """

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
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        smoketest=config["smoketest"],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && python {input.script} harmony {input.data} {params.batch_key} {output} --smoketest={params.smoketest} &> {log}
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
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        smoketest=config["smoketest"],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && python {input.script} pca_harmony {input.data} {params.batch_key} {output} --smoketest={params.smoketest} &> {log}
        """

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
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
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
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && python {input.script} pca_scanorama {input.data} {params.batch_key} {output} &> {log}
        """


rule mnn:
    input:
        data="outputs/{scenario}/{pipeline}.parquet",
        script="scripts/correct_with_mnn.py"
    output:
        "outputs/{scenario}/{pipeline}_mnn.parquet"
    log:
        "logs/{scenario}/{pipeline}_mnn.log"
    conda:
        "../envs/mnn.yaml"  
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && python {input.script} {input.data} {params.batch_key} {output} &> {log}
        """

rule desc:
    input:
        data="outputs/{scenario}/{pipeline}.parquet",
        script="scripts/correct_with_desc.py"
    output:
        "outputs/{scenario}/{pipeline}_desc.parquet"
    log:
        "logs/{scenario}/{pipeline}_desc.log"
    conda:
        "../envs/desc.yaml"  
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && python {input.script} {input.data} {params.batch_key} {output} &> {log}
        """

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
        smoketest=config["smoketest"],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && python {input.script} {input.data} {params.batch_key} {params.label_key} {output} --smoketest={params.smoketest} &> {log}
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
        smoketest=config["smoketest"],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && python {input.script} {input.data} '{params.batch_key}' {params.label_key} {output} --smoketest={params.smoketest} &> {log}
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
        smoketest=config["smoketest"],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && python {input.script} {input.data} '{params.batch_key}' {params.label_key} {output} --smoketest={params.smoketest} &> {log}
        """


rule fastMNN:
    input:
        data="outputs/{scenario}/{pipeline}.parquet",
        script="scripts/correct_with_fastmnn.R"
    output:
        "outputs/{scenario}/{pipeline}_fastMNN.parquet",
    log:
        "logs/{scenario}/{pipeline}_fastmnn.log"
    conda:
        "../envs/fastmnn.yaml"  
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && Rscript {input.script} {input.data} {params.batch_key} {output} &> {log}
        """

rule seurat_cca:
    input:
        data="outputs/{scenario}/{pipeline}.parquet",
        script="scripts/correct_with_seurat.R",
    output:
        "outputs/{scenario}/{pipeline}_seurat_cca.parquet",
    log:
        "logs/{scenario}/{pipeline}_seurat_cca.log"
    conda:
        "../envs/seurat.yaml"  
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && Rscript {input.script} {input.data} {params.batch_key} cca {output} &> {log}
        """


rule seurat_rpca:
    input:
        data="outputs/{scenario}/{pipeline}.parquet",
        script="scripts/correct_with_seurat.R",
    output:
        "outputs/{scenario}/{pipeline}_seurat_rpca.parquet",
    log:
        "logs/{scenario}/{pipeline}_seurat_rpca.log"
    conda:
        "../envs/seurat.yaml"  
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && Rscript {input.script} {input.data} {params.batch_key} rpca {output} &> {log}
        """
