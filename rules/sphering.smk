reg_opts = sphering.log_uniform_sampling(size=config["sphering_n_opts"])


rule sphering_explore:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/sphering/exploration/{pipeline}_reg~{reg}.parquet",
        "outputs/{scenario}/sphering/exploration/{pipeline}_reg~{reg}.npz",
    params:
        mode=config["sphering_mode"],
        reg=lambda wc: float(wc.reg),
        column_norm=config["column_norm"],
        values_norm=config["values_norm"],
    run:
        sphering.sphering(*input, *params, *output)


rule select_best_sphering:
    input:
        parquet_files=expand(
            "outputs/{{scenario}}/sphering/exploration/{{pipeline}}_reg~{reg}.parquet",
            reg=reg_opts,
        ),
        map_files=expand(
            "outputs/{{scenario}}/sphering/exploration/metrics/{{criteria}}/{{pipeline}}_reg~{reg}_map.parquet",
            reg=reg_opts,
        ),
    output:
        parquet_path="outputs/{scenario}/{pipeline}_sphering_{criteria}.parquet",
        map_path="outputs/{scenario}/metrics/{criteria}/{pipeline}_sphering_map.parquet",
    run:
        sphering.select_best(input.parquet_files, input.map_files, *output)


# Avoid mAP recomputation
ruleorder: select_best_sphering > mean_average_precision
