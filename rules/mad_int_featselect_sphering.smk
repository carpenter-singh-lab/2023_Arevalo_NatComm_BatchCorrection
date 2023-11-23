reg_opts = sphering.log_uniform_sampling(size=config["sphering_n_opts"])


rule select_best_map_target2_mad_int_featselect_sphering:
    input:
        map_files=expand(
            "outputs/{{scenario}}/sphering/exploration/map_target2_mad_int_featselect_sphering_reg~{reg}.parquet",
            reg=reg_opts,
        ),
        parquet_files=expand(
            "outputs/{{scenario}}/sphering/exploration/mad_int_featselect_sphering_reg~{reg}.parquet",
            reg=reg_opts,
        ),
    output:
        "outputs/{scenario}/mad_int_featselect_sphering.parquet",
        "outputs/{scenario}/map_target2_mad_int_featselect_sphering.parquet",
    run:
        sphering.select_best(input.map_files, input.parquet_files, *output)


rule mad_int_featselect_sphering_job:
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


rule ap_target2_mad_int_featselect_sphering:
    input:
        "outputs/{scenario}/sphering/exploration/mad_int_featselect_sphering_reg~{reg}.parquet",
    output:
        "outputs/{scenario}/sphering/exploration/ap_target2_mad_int_featselect_sphering_reg~{reg}.parquet",
    params:
        plate_types=["TARGET2"],
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_target2_mad_int_featselect_sphering:
    input:
        "outputs/{scenario}/sphering/exploration/ap_target2_mad_int_featselect_sphering_reg~{reg}.parquet",
    output:
        "outputs/{scenario}/sphering/exploration/map_target2_mad_int_featselect_sphering_reg~{reg}.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)
