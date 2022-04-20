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
def simulated_Kotliarov2020_processed(n_var = 10, n_obs = 100):
    #adatam = read('./data/Kotliarov2020_processed_citeseq_merged_annotated.h5ad', cache=False)
    #print(adata.obs["leiden"])
    raw_data = [[y if bool(random.getrandbits(1)) else -abs(y) for y in x] for x in np.random.exponential(1.6, size=(n_obs, n_var))]
    counts = csr_matrix(raw_data, dtype=np.float32)
    adata = AnnData(counts)
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
    """ count_variable = 'leiden'
    subset_variable = 'batch'

    subsets = natsorted(adata.obs.get(subset_variable).value_counts().index.tolist())

    for subset in subsets:
        datam = _subset_adata(
            adatam, filter_criteria=adatam.obs.get(subset_variable) == subset, raw=False
        )
        all_countsm = count_occurrence(
            datam, count_variable=count_variable, add_percentage=True
        )
        data = _subset_adata(
            adata, filter_criteria=adata.obs.get(subset_variable) == subset, raw=False
        )
        all_counts = count_occurrence(
            data, count_variable=count_variable, add_percentage=True
        )
        print(all_counts)
        print(all_countsm) """
    #bc.tl.count_occurrence_subset(
    #    adata, subset_variable, count_variable=count_variable, return_percentage=False
    #)
    #print(adata.obs["leiden"])
    #print(adata)
    return adata
from pandas import DataFrame, melt
""" def count_occurrence(adata, count_variable="celltype", add_percentage=False):
    data = adata.obs.get(count_variable)
    if data is None:
        sys.exit(
            "please specify a column name (count_variable) that is present in adata.obs"
        )
    print(data)
    print('super duper')
    # get counts for specified column
    counts = adata.obs.get(count_variable).value_counts()
    print(type(counts))
    total_counts = sum(counts.tolist())
    print(total_counts)
    # initialize a dataframe to store information in
    data = DataFrame(index=counts.index.tolist(), data={"Counts": counts.tolist()})
    print('dataframe has been created')
    if add_percentage:  ##  calculate percentages
        percentages = []
        print(data)
        print(counts)
        for i in range(len(counts)):
            print(i)
            count = data["Counts"][i]
            percentages.append(round(count / total_counts * 100, 2))
        # add percentages to dataframe
        data["Percentage"] = percentages
    print('returning data')
    return data
 """
"python -c 'from _helper import *; print(simulated_pbmc3k_raw())'"