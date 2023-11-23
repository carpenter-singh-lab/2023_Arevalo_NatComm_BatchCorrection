reg_opts = sphering.log_uniform_sampling(size=config["sphering_n_opts"])


rule select_best_map_mad_int_featselect_sphering:
    input:
        map_files=expand(
            "outputs/{{scenario}}/sphering/exploration/map_{{criteria}}_mad_int_featselect_sphering_reg~{reg}.parquet",
            reg=reg_opts,
        ),
        parquet_files=expand(
            "outputs/{{scenario}}/sphering/exploration/mad_int_featselect_sphering_reg~{reg}.parquet",
            reg=reg_opts,
        ),
    output:
        map_path="outputs/{scenario}/map_{criteria}_mad_int_featselect_sphering.parquet",
    params:
        # Required to be in params due to the lack of criteria wildcard
        parquet_path=lambda wc: f"outputs/{wc.scenario}/mad_int_featselect_sphering.parquet",
    run:
        sphering.select_best(input.map_files,
                             input.parquet_files,
                             params.parquet_path,
                             output.map_path)

# Avoid mAP recomputation
ruleorder: select_best_map_mad_int_featselect_sphering > mean_average_precision

rule mad_int_featselect_sphering_explore:
    input:
        "outputs/{scenario}/mad_int_featselect.parquet",
    output:
        "outputs/{scenario}/sphering/exploration/mad_int_featselect_sphering_reg~{reg}.parquet",
        "outputs/{scenario}/sphering/exploration/mad_int_featselect_spherer_reg~{reg}.npz",
    params:
        mode=config["sphering_mode"],
        reg=lambda wc: float(wc.reg),
        column_norm=config["column_norm"],
        values_norm=config["values_norm"],
    run:
        sphering.sphering(*input, *params, *output)
