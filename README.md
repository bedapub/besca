# BESCA (BEDA's Single Cell Analysis library)

The BESCA (BEDA’s single cell sequencing analysis) package contains many useful python functions to use for your single-cell analysis.

The package has been grouped into 3 categories:  

- preprocessing functions: this submodule contains all functions relevant to data preprocessing  
- plotting functions: additional plot types not available in the standard scanpy package  
- tools: contains additional tools to e.g. perform differential gene analysis or load/export data  

For more information please view the package documentation: https://bedapub.github.io/besca/



## Installation

If you are familiar with python packages simply install it using pip:  

```
pip install git+https://github.com/bedapub/besca.git
```

Besca comes with a binary called reformat written in C and was compiled in linux-64. Therefore, besca runs exclusively on linux-64.
In some cases the binary needs to be made executeable. To do so show the location of the path and navigate to the besca package.  
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

If you are not very familiar with python packages here is a detailled description.  

If you don't have a conda python installation download and install [miniconda](https://docs.conda.io/en/latest/miniconda.html). While installing I recommend to accept everything asked by the miniconda installation.  

As a next step we create a separate environment for besca which is also called besca.  
```
conda create --name besca python=3.7.1
```  
We can activate this environment.  
```
conda activate besca
```
Within this enviroment we are gonna install besca using pip.  
```
pip install git+https://github.com/bedapub/besca.git
```
We need to make a binary which comes with besca executeable. To do so show the location of the path and navigate to the besca package.  
```
pip show besca
cd Location/besca
```
Navigate in the directory with the binary and make it executeable.  
```
cd export
chmod u+x reformat
```

You should now have successfully installed besca.

## Running besca on a HPC with a SLURM worklod manager  

If you have access to an HCP which uses SLURM as a workload manger you can run the jupyter notebooks coming with besca located in `workbooks/` with dedicated resources.  
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