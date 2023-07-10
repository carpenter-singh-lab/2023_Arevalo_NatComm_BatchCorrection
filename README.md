# Batch correction baseline

Scripts to evaluate batch correction methods in Target2 plates. This script generates several embedding visualizations and .csv files for quantitative evaluation of batch correction methods.

## Installation

This set of commands creates the environment from scratch and install the packages via `pip` command.

```bash
./Miniconda3-py39_4.11.0-Linux-x86_64.sh
conda update -n base -c defaults conda
conda install mamba -n base -c conda-forge
mamba init
mamba create -n bcorrect python=3.10 R=3.6 -c conda-forge
mamba activate bcorrect
mamba install r-devtools
mamba install numpy pandas scvi-tools tensorflow-probability pyarrow plotly scikit-learn sqlalchemy  pytest leidenalg umap-learn numba anndata  scanpy scikit-misc deprecated  llvmlite louvain pillow pydot kneed python-igraph rpy2 pymde python-kaleido
mamba install anndata2ri -c bioconda

# Install correction methods
pip install harmonypy pymde scanorama desc
# pip install mnnpy. # requires installation from source
```

The following packages are installed manually

```bash
git clone https://github.com/cytomining/pycytominer ~/projects/pycytominer
cd ~/projects/pycytominer
pip install -e .

git clone git@github.com:theislab/scib.git ~/projects/scib
cd ~/projects/scib
pip install -e .

pip install git+https://github.com/cytomining/copairs.git@v0.0.3
```

Additionally, in order to run the R package `kBET`, you need to install it through R.

```R
devtools::install_github('theislab/kBET')
```


## Get input data

Download profiles and metadata:
```bash
source download_data.sh
```

## Setup file

The [`config.json`](config.json) file has the parameters to define the sources to be the processed as well as the preprocessing steps:

 - `"sources": ["source_4"]`: List of sources to be processed.
 - `"plate_types": ["TARGET2"]`: plate types to be loaded.
 - `"epsilon_mad": 1e0`: Hyperparameter used in the Robust MAD Normalization. Only used when `mad_norm` set to `true`.
 - `"mad_norm": true`: Enable the Robust MAD normalization.
 - `"outlier_removal": true`: Enable the outlier removal preprocessing.
 - `"nan_removal": true`: Enable the NaN removal preprocessing.
 - `"feature_selection": true`: Enable the feature selection preprocessing.
 - `"sphering": true`: Enable the sphering alignment.
 - `"sphering_mode": "corr"`: Hyperparameter used to define how the normalization matrix is calculated in the sphering.
 - `"sphering_lambda": null`: Hyperparameter to regularize sphering when it is set to `null` the parameter is auto configured based on the eigenvalues.

## Run the analysis
To train all the available models, run:
```bash
python run_correction config.json
```

Output is written in the `./outputs/` folder.
