"""Simulated datasets.
"""

from anndata import AnnData
import numpy as np
import besca as bc
from scipy.sparse import csr_matrix
from string import ascii_uppercase
import random

def simulated_pbmc3k_raw(n_var = 100, n_obs = 10):
    counts = csr_matrix(np.random.poisson(1, size=(n_obs, n_var)), dtype=np.float32)
    adata = AnnData(counts)
    obs_names = [''.join((random.choice('ACGT') for x in range(14))) for i in range(adata.n_obs)]
    obs_names = [''.join(f"{name}-1") for name in obs_names]
    adata.obs_names = obs_names
    adata.obs["CELL"] = obs_names
    adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
    adata.var["SYMBOL"] = [f"Gene_{i:d}" for i in range(adata.n_vars)]
    return adata


import pkg_resources
from scanpy import read
from natsort import natsorted
from besca._helper import subset_adata as _subset_adata
from pandas.api.types import CategoricalDtype
def simulated_Kotliarov2020_processed(n_var = 25, n_obs = 250):
    #adatam = read('./data/Kotliarov2020_processed_citeseq_merged_annotated.h5ad', cache=False)
    #adatam.X = adatam.X.toarray()
    #iLoc = adatam.var_names.tolist().index('CD19')
    #print(adatam.X[iLoc].flatten())
    raw_data = [[y if bool(random.getrandbits(1)) else -abs(y) for y in x] for x in np.random.exponential(1.6, size=(n_obs, n_var))]
    counts = csr_matrix(raw_data, dtype=np.float32)
    raw_adata = AnnData(counts)
    raw_adata.var_names = [f"Gene_{i:d}" for i in range(raw_adata.n_vars)]
    adata = AnnData(counts)
    adata.raw = raw_adata
    adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
    obs_names = [''.join((random.choice('ACGT') for x in range(14))) for i in range(adata.n_obs)]
    obs_names = [''.join(f"10X_CiteSeq_donor{random.randint(1,10)}.{name}-1") for name in obs_names]
    adata.obs_names = obs_names
    adata.obs["CELL"] = obs_names
    adata.obs["batch"] =  np.random.choice([1,2], adata.n_obs)
    adata.obs["leiden"] =  [str(random.randint(0,22)) for x in range(adata.n_obs)]
    status_type = CategoricalDtype(categories=[''.join(f"{x}") for x in range(23)], ordered=True)
    adata.obs["leiden"] = adata.obs["leiden"].astype(status_type)
    adata.obs["donor"] =  [''.join(f"donor{random.randint(1,10)}") for x in range(adata.n_obs)]
    adata.obs["donor"] = adata.obs["donor"].astype("category")
    adata.obs["CONDITION"] = ["PBMC_healthy" for x in range(adata.n_obs)]
    adata.obs["celltype0"] = ["hematopoietic cell" for x in range(adata.n_obs)]
    adata.obs["sampleid"] = [x for x in range(adata.n_obs)]
    return adata

"python -c 'from _helper import *; print(simulated_pbmc3k_raw())'"