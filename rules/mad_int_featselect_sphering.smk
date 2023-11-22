reg_opts = correct.log_uniform_sampling(size=config["sphering_n_opts"])


rule select_best:
    input:
        map_files=expand(
            "outputs/{{scenario}}/sphering/exploration/reg~{reg}/map_target2_mad_int_featselect_sphered.parquet",
            reg=reg_opts,
        ),
        parquet_files=expand(
            "outputs/{{scenario}}/sphering/exploration/reg~{reg}/mad_int_featselect_sphered.parquet",
            reg=reg_opts,
        ),
    output:
        "outputs/{scenario}/sphering/mad_int_featselect_sphered.parquet",
    run:
        correct.select_best(input.map_files, input.parquet_files, *output)


rule mad_int_featselect_sphering_job:
    input:
        "outputs/{scenario}/mad_int_featselect.parquet",
    output:
        "outputs/{scenario}/sphering/exploration/reg~{reg}/mad_int_featselect_sphered.parquet",
        "outputs/{scenario}/sphering/exploration/reg~{reg}/mad_int_featselect_spherer.npz",
    params:
        mode=config["sphering_mode"],
        reg=lambda wc: float(wc.reg),
        column_norm=config["column_norm"],
        values_norm=config["values_norm"],
    run:
        correct.sphering(*input, *params, *output)


rule ap_target2_mad_int_featselect_sphering:
    input:
        "outputs/{scenario}/sphering/exploration/reg~{reg}/mad_int_featselect_sphered.parquet",
    output:
        "outputs/{scenario}/sphering/exploration/reg~{reg}/ap_target2_mad_int_featselect_sphered.parquet",
    params:
        plate_types=["TARGET2"],
    run:
        qc.metrics.average_precision(*input, *output, **params)


rule map_target2_mad_int_featselect_sphering:
    input:
        "outputs/{scenario}/sphering/exploration/reg~{reg}/ap_target2_mad_int_featselect_sphered.parquet",
    output:
        "outputs/{scenario}/sphering/exploration/reg~{reg}/map_target2_mad_int_featselect_sphered.parquet",
    run:
        qc.metrics.mean_average_precision(*input, *output)
