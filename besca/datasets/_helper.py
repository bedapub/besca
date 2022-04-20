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

"python -c 'from _helper import *; print(simulated_pbmc3k_raw())'"