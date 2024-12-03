configfile: "./inputs/conf/scenario_6.json"


# Define processing workflows and correction methods to run
WORKFLOWS = [
    "mad_int_featselect",
]

METHODS = [
    #"scanorama", problem
    #"fastMNN", problem latency
    #"mnn", problem latency
    #"harmony", 
    #"combat", 
    #"desc", 
    #"scvi", 
    #"scanvi",
    "gaushvi",
    "gaushanvi", 
    # "sysvi", problem
    #"scpoli", problem
    #"sphering", 
    #"seurat_cca", problem latency
    #"seurat_rpca", problem latency
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
