# rule average_precision_negcon:
#     input:
#         "outputs/{scenario}/{pipeline}.parquet",
#     output:
#         "outputs/{scenario}/metrics/{criteria}/{pipeline}_ap_negcon.parquet",
#     params:
#         plate_types=lambda wc: ["TARGET2"] if wc.criteria == "target2" else ["COMPOUND"],
#     run:
#         metrics.average_precision_negcon(*input, *output, **params)


# rule average_precision_nonrep:
#     input:
#         "outputs/{scenario}/{pipeline}.parquet",
#     output:
#         "outputs/{scenario}/metrics/{criteria}/{pipeline}_ap_nonrep.parquet",
#     params:
#         plate_types=lambda wc: ["TARGET2"] if wc.criteria == "target2" else ["COMPOUND"],
#     run:
#         metrics.average_precision_nonrep(*input, *output, **params)


