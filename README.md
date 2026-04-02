# BESCA (BEyond Single Cell Analysis)

[![Run doctests](https://github.com/bedapub/besca/actions/workflows/doc-tests.yml/badge.svg)](https://github.com/bedapub/besca/actions/workflows/doc-tests.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

The BESCA (BEyond Single Cell Analysis) package contains many useful python functions to use for your single-cell analysis.

The package has been grouped into 3 categories:  

- preprocessing functions: this submodule contains all functions relevant to data preprocessing  
- plotting functions: additional plot types not available in the standard scanpy package  
- tools: contains additional tools to e.g. perform differential gene analysis or load/export data  

For more information please view the package documentation: https://bedapub.github.io/besca/

Please consider citing our publication if you use Besca for your research:

- Mädler SC, Julien-Laferriere A, Wyss L, Phan M, Sonrel A, Kang ASW, Ulrich E, Schmucki R, Zhang JD, Ebeling M, Badi L, Kam-Thong T, Schwalie PC, Hatje K. <a href="https://doi.org/10.1093/nargab/lqab102" target="_blank">Besca, a single-cell transcriptomics analysis toolkit to accelerate translational research</a>. <i>NAR Genom Bioinform</i>. 2021


If you are interested in contributing you can check the repository wiki for helpful information on contributing: https://github.com/bedapub/besca/wiki

For faster/smaller download, in case of slow internet connection or low storage capacity, please use following command to clone this repository:
```
git clone --filter=blob:none git@github.com:bedapub/besca.git
```

## Installation

From version 3.0.0, Besca requires Python 3.12 or above.

If you are familiar with python packages simply install them using pip: 

```
pip install besca
```

or

```
pip install git+https://github.com/bedapub/besca.git
```

### Python beginner guide

If you are not very familiar with python packages here is a detailed description.  

If you don't have a conda python installation download and install [miniconda](https://docs.conda.io/en/latest/miniconda.html). While installing we recommend accepting everything asked by the miniconda installation.  

As a next step, we create a separate environment for besca which is also called besca.  

```
conda create --name besca python=3.12
```  

We can activate this environment.  

```
conda activate besca
```

Within this environment, we can install besca using pip.  

```
pip install git+https://github.com/bedapub/besca.git
```

You should now have successfully installed besca.

In case you met any problems, please report an issue.

To install [Jupyter Notebook](https://jupyter.readthedocs.io/en/latest/install/notebook-classic.html), type

```
conda install jupyter
```

and type 

```
jupyter notebook
```

to start a Jupyter Notebook in your browser. See [documentation](https://jupyter.readthedocs.io/en/latest/running.html#running) for further details. 


## Running besca on an HPC with a SLURM workload manager  

If you have access to an HPC which uses SLURM as a workload manager you can run the jupyter notebooks coming with besca located in `workbooks/` with dedicated resources.  
To do so, start an interactive session on your HPC.  

```
interactive -c 8 -m 16G -t 180 # This allocates 8 CPUs, 16 GB of memory for 3 hours
```

If you have installed besca in a conda environment like explained above activate the environment.  

```
conda activate besca
```

Start a jupyter notebook.  

```
jupyter-notebook --ip=* --no-browser
```

You can now run the jupyter notebooks coming with besca.



## Datasets and Analysis notebooks


Besca run-examples and datasets annotation notebooks can be found in:


[https://github.com/bedapub/besca_publication_results](https://github.com/bedapub/besca_publication_results)


All processed datasets were uploaded to Zenodo, within the Besca community:

[https://zenodo.org/communities/besca/](https://zenodo.org/communities/besca/)






