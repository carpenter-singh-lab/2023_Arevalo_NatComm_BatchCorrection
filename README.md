# Evaluating batch correction methods for image-based cell profiling

This repository contains the source code to reproduce the results in the
paper: "Evaluating batch correction methods for image-based cell profiling".

These scripts generate several embedding visualizations and .csv files for
quantitative evaluation of batch correction methods.

## Installation

We suggest [Mamba](https://github.com/conda-forge/miniforge#mambaforge) for
environment management. The following commands create the environment from
scratch and install the required packages.

```bash
mamba env create --file environment.yaml
mamba activate batchcp
```
### kBET installation
Run the following command to install R package `kBET`:

```bash
R -e "remotes::install_github('theislab/kBET')"
```

### scib installation

The single-cell integration benchmark `scib` package includes `c/c++` code that
should be compiled specifically for your environment. Thus it requires
installation from source. More info in [this
issue](https://github.com/theislab/scib/issues/308)

```bash
DEST=$HOME/projects/scib
git clone https://github.com/theislab/scib.git $DEST
cd $DEST
git checkout v1.1.4
pip install -e .
```

### mnnpy installation
Similarly, `mnnpy` may require manual installation. More info at
https://github.com/chriscainx/mnnpy#install

## Get input data

Download profiles and metadata:
```bash
source download_data.sh
```

## Run scenarios
Every scenario reported in the paper can be reproduced running the associated
bash script. For example, to reproduce Scenario 1:

```bash
source run_scenario_1.sh
```

The output is written in the `./outputs/` folder. You can get the scores and
plots in [Plot_figures.ipynb](Plot_figures.ipynb) notebook.
