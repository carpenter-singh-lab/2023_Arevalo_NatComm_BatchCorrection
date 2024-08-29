reg_opts = sphering.log_uniform_sampling(size=25)


rule sphering_explore:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/sphering/exploration/{pipeline}_reg~{reg}.parquet",
        "outputs/{scenario}/sphering/exploration/{pipeline}_reg~{reg}.npz",
    params:
        method=config["sphering_method"],
        reg=lambda wc: float(wc.reg),
        column_norm=config["column_norm"],
        values_norm=config["values_norm"],
    run:
        sphering.sphering(*input, *params, *output)


rule select_best_sphering:
    input:
        parquet_files=expand(
            f"outputs/{scenario}/sphering/exploration/{{pipeline}}_reg~{{reg}}.parquet",
            reg=reg_opts,
            allow_missing=True,
        ),
        map_negcon_files=expand(
            f"outputs/{scenario}/sphering/exploration/metrics/{criteria}/{{pipeline}}_reg~{{reg}}_map_negcon.parquet",
            reg=reg_opts,
            allow_missing=True,
        ),
        map_nonrep_files=expand(
            f"outputs/{scenario}/sphering/exploration/metrics/{criteria}/{{pipeline}}_reg~{{reg}}_map_nonrep.parquet",
            reg=reg_opts,
            allow_missing=True,
        ),
    output:
        parquet_path=f"outputs/{scenario}/{{pipeline}}_sphering.parquet",
        ap_negcon_path=f"outputs/{scenario}/metrics/{criteria}/{{pipeline}}_sphering_ap_negcon.parquet",
        ap_nonrep_path=f"outputs/{scenario}/metrics/{criteria}/{{pipeline}}_sphering_ap_nonrep.parquet",
        map_negcon_path=f"outputs/{scenario}/metrics/{criteria}/{{pipeline}}_sphering_map_negcon.parquet",
        map_nonrep_path=f"outputs/{scenario}/metrics/{criteria}/{{pipeline}}_sphering_map_nonrep.parquet",
    run:
        sphering.select_best(
            input.parquet_files,
            input.map_negcon_files,
            input.map_nonrep_files,
            output.ap_negcon_path,
            output.ap_nonrep_path,
            output.map_negcon_path,
            output.map_nonrep_path,
            output.parquet_path,
        )


# Because map files
ruleorder: select_best_sphering > mean_average_precision
ruleorder: select_best_sphering > average_precision_nonrep
ruleorder: select_best_sphering > average_precision_negcon