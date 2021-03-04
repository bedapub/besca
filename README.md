# BESCA (BEDA's Single Cell Analysis library)

The BESCA (BEDA’s single cell sequencing analysis) package contains many useful python functions to use for your single-cell analysis.

The package has been grouped into 3 categories:  

- preprocessing functions: this submodule contains all functions relevant to data preprocessing  
- plotting functions: additional plot types not available in the standard scanpy package  
- tools: contains additional tools to e.g. perform differential gene analysis or load/export data  

For more information please view the package documentation: https://bedapub.github.io/besca/

Please find our preprint posted on bioRxiv here: https://biorxiv.org/cgi/content/short/2020.08.11.245795v2

If you are interested in contributing you can check the repository wiki for helpful information on contributing: https://github.com/bedapub/besca/wiki

## Installation

If you are familiar with python packages simply install them using pip:  

```
pip install git+https://github.com/bedapub/besca.git
```

Besca comes with a binary called reformat written in C and was compiled in linux-64. Therefore, besca runs exclusively on linux-64.


### Set the executable flag to the binary file `reformat` <a name="binary"></a>

In some cases, the binary file needs to be made executable. To do so, run the following one-liner.

```bash
pip show besca | grep Location | cut -f 2 -d ":" | awk -v OFS="" '{print "chmod u+x" $0 "/besca/export/reformat"}' | bash
```

If you want to avoid piping to bash, or want to it step by step, here is how to. Show the location of the path and navigate to the besca package.  

```
pip show besca
cd Location/besca
```

Navigate in the directory containing the binary and make it executeable.  

```
cd export
chmod u+x reformat
```

### Python beginner guide

If you are not very familiar with python packages here is a detailed description.  

If you don't have a conda python installation download and install [miniconda](https://docs.conda.io/en/latest/miniconda.html). While installing we recommend accepting everything asked by the miniconda installation.  

As a next step, we create a separate environment for besca which is also called besca.  

```
conda create --name besca python=3.7.1
```  

We can activate this environment.  

```
conda activate besca
```

Within this environment, we can install besca using pip.  

```
pip install git+https://github.com/bedapub/besca.git
```

Now following the [instruction above](#binary) to set the executable flag to the binary file shipped with besca.

You should now have successfully installed besca.

In case you met any problems, please report an issue.


### R dependencies for additional methods

Although the standard workflow can be run without any R dependencies, BESCA can run a selection of performant methods developed in R. These additional methods are :

- [`isOutlier`](https://www.rdocumentation.org/packages/scater/versions/1.0.4/topics/isOutlier) from `scater`: for outlier detection and filtering recommendations. Implemented in the `besca.pp.valOutlier` function.  
- [`SCTransform`](https://rdrr.io/github/satijalab/seurat/man/SCTransform.html) : one of the normalization methods proposed by the `Seurat` package. Implemented in the `besca.pp.scTransform` function. 
- [`maxLikGlobalDimEst`](https://cran.r-project.org/web/packages/intrinsicDimension/intrinsicDimension.pdf) from `intrinsicDimension` : for an estimation of the number of dimensions to use for clustering. Implemented in the `besca.st.maxLikGlobalDimEst` function. 
- [`deviance`](https://rdrr.io/bioc/scry/man/devianceFeatureSelection.html) and [`VST`](https://rdrr.io/github/satijalab/seurat/man/SCTransform.html): for highly-variable genes selection. Implemented in the `besca.st.deviance` function. 
- [`DSB`](https://github.com/niaid/dsb): for denoising ADT counts data based on background noise. Implemented in the `besca.st.dsb_normalize` function.  


#### Conda installation

If you used a conda enviroment it is possible to install most needed dependencies using Conda too. 

With an activated environment using:

```
conda activate besca
```


One can run the commands below:

```

conda install -y -c conda-forge r=4.0 rpy2 r-essentials r-base r-devtools r-withr r-vctrs r-tidyverse r-magrittr r-data.table r-Matrix r-ggplot2 r-readr r-seurat r-intrinsicdimension r-mclust r-sitmo r-patchwork --force-install
conda install -y -c bioconda anndata2ri R bioconductor-dropletutils bioconductor-scry
conda install -c bioconda bioconductor-scater
```

This should install in your conda envrionment the dependencies under : *conda_path/lib/R/library* of your conda environment path.

#### Pip installation


If you want to run one of these methods in the workflow, please install the required libraries by running the following commands in the `besca` installation directory (or simply download the `Rlibs.R` file):

```
pip install rpy2 anndata2ri
<conda_bin_path>/Rscript Rlibs.R <conda_R_library_path>
 ```
### Location of conda R bin
`~/.conda/envs/[environnement_name]/bin/`

### Location of the conda R library 

If you used conda, by default, libraries should be installed into your conda environment path, typically `~/.conda/envs/[environnement_name]/lib/R/library`.
If this is not the right path, please verify the path to your conda enviroment using `conda list env`.


To minimize risks conflicts between libraries, it is advised to set your `your_R_library_path` to such path also while using pip.
 
In the standard workflow notebook, all of these methods are controlled through the `r_methods` option but it is of course possible to manually switch between them and the standard workflow. Please also specify the location of your R library with the `rlib_loc` option of the notebook.  

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






