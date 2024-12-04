# Evaluating batch correction methods for image-based cell profiling

This repository adapts the source code from the paper: "Evaluating batch correction methods for image-based cell profiling".


## Installation

We suggest [Mamba](https://github.com/conda-forge/miniforge#mambaforge) for
environment management. The following commands create the environment from
scratch and install the required packages.

```bash
mamba env create --file environment.yaml
mamba activate batchcp
```

or if you already have an environment you just want to update:

```bash
mamba env update --name batchcp --file environment.yaml --prune
```

### mnnpy installation

Similarly, `mnnpy` may require manual installation. More info at
https://github.com/chriscainx/mnnpy#install

## Get input data

Download profiles and metadata:
```bash
bash download_data.sh
```

## Run scenarios
Every scenario reported in the paper can be reproduced running snakemake with
the associated config file. For example, to reproduce Scenario 6 using 8 cores:

```bash
snakemake -c8 --configfile inputs/conf/scenario_7.json --use-conda --conda-prefix "./env_store/" --resources nvidia_gpu=2
```

```bash
snakemake -c8 --configfile inputs/conf/scenario_7.json --use-conda --conda-prefix "./env_store/" --rulegraph > rulegraph.dot
```

snakemake -c8 --configfile inputs/conf/scenario_6.json --use-conda --conda-prefix "./env_store/" --until results_table

You can get the scores, corrected profiles and plots in the `./outputs` folder.


## Learnings
- constraining GPU ressources is not worth it because we're not maxing out the VRAM for most methods on target2