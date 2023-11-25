rule harmony:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_harmony.parquet",
    params:
        batch_key=config['batch_key']
    run:
        correct.harmony(*input, *params, *output)
