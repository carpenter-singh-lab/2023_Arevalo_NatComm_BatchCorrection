rule get_labels:
    params:
        drugs_url="https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt",
        inchikey_mapper_url="https://raw.githubusercontent.com/jump-cellpainting/JUMP-Target/fc1561685b10d98543b8ba5d2484e960778324a9/JUMP-Target-2_compound_metadata.tsv",
        compound_path="inputs/metadata/compound.csv.gz",
    output:
        "inputs/metadata/target_annotations.parquet",
    run:
        metrics.consistency.get_labels(*params, *output)


rule median_consensus:
    input:
        "outputs/{prefix}/{pipeline}.parquet",
    output:
        "outputs/{prefix}/median_consensus/{pipeline}.parquet",
    run:
        metrics.consistency.median_profile(*input, *output)


rule annotate_consensus:
    input:
        "outputs/{prefix}/median_consensus/{pipeline}.parquet",
        "inputs/metadata/target_annotations.parquet",
    output:
        "outputs/{prefix}/metrics/target2/consistency/profiles/{pipeline}.parquet",
    run:
        metrics.consistency.annotate_median_profile(*input, *output)


rule target_ap:
    input:
        "outputs/{prefix}/metrics/target2/consistency/profiles/{pipeline}.parquet",
    output:
        "outputs/{prefix}/metrics/target2/consistency/{pipeline}_ap.parquet",
    run:
        metrics.consistency.target_ap(*input, *output)


rule target_map:
    input:
        "outputs/{prefix}/metrics/target2/consistency/{pipeline}_ap.parquet",
    output:
        "outputs/{prefix}/metrics/target2/consistency/{pipeline}_map.parquet",
    run:
        metrics.consistency.target_map(*input, *output)
