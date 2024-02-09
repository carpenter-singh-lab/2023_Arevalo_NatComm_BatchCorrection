configfile: "./inputs/conf/scenario_1.json"


# Define processing workflows and correction methods to run
WORKFLOWS = [
    "mad_int_featselect",
]
METHODS = ["scanorama", "mnn", "harmony", "combat", "desc", "scvi", "sphering"]


# Load rules
include: "rules/common.smk"
include: "rules/processing.smk"
include: "rules/metrics.smk"
include: "rules/correct.smk"
include: "rules/projection.smk"
include: "rules/plots.smk"


rule all:
    input:
        expand(plots_pattern, plot=PNG_PLOTS, ext=["png"]),
        expand(plots_pattern, plot=SVG_PLOTS, ext=["svg"]),
