
rule optimize_scpoli:
    input:
        data="outputs/{scenario}/" + config["preproc"] + ".parquet",
        script="scripts/optimise_scpoli.py"
    output:
        path="outputs/{scenario}/optimization/optuna_scpoli.csv"
    log:
        "logs/{scenario}/" + config["preproc"] + "_optimize_scpoli.log"
    conda:
        "../envs/scpoli.yaml"
    params:
        batch_key=config["batch_key"] if isinstance(config["batch_key"], str) else config["batch_key"][0],
        label_key=config["label_key"],
        trials=config["optuna_trials"],
        smoketest="--smoketest" if config["smoketest"] else "",
    resources:
        nvidia_gpu=1
    shell:
        """
        export PYTHONPATH=$(dirname $(pwd)):$(pwd) && \
        python '{input.script}' \
            --input_data '{input.data}' \
            --batch_key '{params.batch_key}' \
            --label_key '{params.label_key}' \
            --n_trials '{params.trials}' \
            --output_path '{output.path}' \
            {params.smoketest} \
            &> '{log}'
        """
