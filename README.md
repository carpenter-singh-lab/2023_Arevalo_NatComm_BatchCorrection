# Evaluating batch correction methods for image-based cell profiling

This repository adapts the source code from the paper: "Evaluating batch correction methods for image-based cell profiling".


## Installation

We suggest [Mamba](https://github.com/conda-forge/miniforge#mambaforge) for
environment management. The following commands create the environment from
scratch and install the required packages.

```bash
mamba env create --file environment.yaml
mamba activate batchcp
mamba install pyarrow -c conda-forge
```

or if you already have an environment you just want to update:

```bash
mamba env update --name batchcp --file environment.yaml --prune
```

### kBET installation
Run the following command to install R package `kBET`:

```bash
R -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes'); remotes::install_github('theislab/kBET')"
```

### scib troubleshooting

The single-cell integration benchmark `scib` package may fail because it
includes `c/c++` code that should be compiled specifically for your
environment. An alternative is to install it from source:

```bash
DEST=$HOME/projects/scib
git clone https://github.com/theislab/scib.git $DEST
cd $DEST
git checkout v1.1.4
pip install -e .
```
More info in [this issue](https://github.com/theislab/scib/issues/308)

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
snakemake -c8 --configfile inputs/conf/scenario_6.json --use-conda --conda-prefix "./env_store/"
```

snakemake -c8 --configfile inputs/conf/scenario_6.json --use-conda --conda-prefix "./env_store/" --until tidy_scores

You can get the scores, corrected profiles and plots in the `./outputs` folder.
