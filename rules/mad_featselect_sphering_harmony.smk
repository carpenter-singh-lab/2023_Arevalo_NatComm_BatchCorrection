rule mad_featselect_sphering_harmony:
    input:
        "outputs/{scenario}/mad_featselect_sphering.parquet",
    output:
        "outputs/{scenario}/mad_featselect_sphering_harmony.parquet",
    params:
        batch_key=config['batch_key']
    run:
        correct.harmony(*input, *params, *output)
