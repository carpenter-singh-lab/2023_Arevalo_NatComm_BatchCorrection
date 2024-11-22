configfile: "./inputs/conf/scenario_6.json"


# Define processing workflows and correction methods to run
WORKFLOWS = [
    "mad_int_featselect",
]
METHODS = [
    "scanorama",
    "scanorama_pca",
    # "mnn", # pain in the ass, we'll just use the R version fastMNN instead
    "fastMNN",
    "harmony",
    # "harmony_pca", # performs the same as normal harmony
    "combat",
    "desc",
    "scvi",
    "sysvi",
    "scpoli",
    # "scpoli_pca", # performs quite a bit worse than normal scpoli
    "sphering",
    "seurat_cca",
    "seurat_rpca",
]


# Load rules
include: "rules/common.smk"
include: "rules/processing.smk"
include: "rules/metrics.smk"
include: "rules/correct.smk"
include: "rules/projection.smk"
include: "rules/plots.smk"


rule all:
    input:
        expand(plots_pattern, plot=PDF_PLOTS, ext=["pdf"]),
