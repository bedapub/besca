"""Simulated datasets.
"""

from anndata import AnnData
import numpy as np
import besca as bc
from scipy.sparse import csr_matrix
import random
from pandas.api.types import CategoricalDtype

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




def simulated_Kotliarov2020_processed(n_var = 25, n_obs = 250):
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

def simulated_pbmc3k_processed(n_var = 25, n_obs = 250):
    raw_data = [[y if bool(random.getrandbits(1)) else -abs(y) for y in x] for x in np.random.exponential(1.6, size=(n_obs, n_var))]
    counts = csr_matrix(raw_data, dtype=np.float32)
    raw_adata = AnnData(counts)
    obs_names = [''.join((random.choice('ACGT') for x in range(14))) for i in range(raw_adata.n_obs)]
    obs_names = [''.join(f"{name}-1") for name in obs_names]
    raw_adata.obs_names = obs_names
    raw_adata.obs["CELL"] = obs_names
    raw_adata.var_names = [f"Gene_{i:d}" for i in range(raw_adata.n_vars)]
    raw_adata.var["SYMBOL"] = [f"Gene_{i:d}" for i in range(raw_adata.n_vars)]
    adata = AnnData(np.asarray(raw_data))
    adata.raw = raw_adata
    adata.obs_names = obs_names
    adata.obs["CELL"] = obs_names
    adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
    adata.var["SYMBOL"] = [f"Gene_{i:d}" for i in range(adata.n_vars)]
    adata.obs["dblabel"] = np.random.choice(["naive B cell", "central memory CD4-positive, alpha-beta T cell", "classical monocyte"], size=(adata.n_obs,))
    adata.obs["leiden"] =  [str(random.randint(0,12)) for x in range(adata.n_obs)]
    status_type = CategoricalDtype(categories=[f"{x}" for x in range(13)], ordered=True)
    adata.obs["leiden"] = adata.obs["leiden"].astype(status_type)
    return adata
""" import scanpy as sc
from besca._helper import subset_adata as _subset_adata
from scanpy.preprocessing import highly_variable_genes as sc_highly_variable_genes
import logging
from scanpy.preprocessing import scale as sc_scale
from scanpy.tools import pca as sc_pca
from scanpy.preprocessing import filter_genes_dispersion, log1p, regress_out, neighbors
from scanpy.tools import umap, louvain, leiden
def recluster(
    adata,
    celltype,
    celltype_label="leiden",
    min_mean=0.0125,
    max_mean=4,
    min_disp=0.5,
    resolution=1.0,
    regress_out_key=None,
    random_seed=0,
    show_plot_filter=False,
    method="leiden",
    batch_key=None,
    n_shared=2,
):

    >>> import besca as bc
    >>> import scanpy as sc
    >>> adata = bc.datasets.simulated_pbmc3k_processed()
    >>> adata_subset = bc.tl.rc.recluster(adata, celltype=('0', '1', '3', '6'), resolution = 1.3)
    >>> sc.pl.umap(adata_subset, color = ['leiden', 'Gene_4', 'Gene_5', 'Gene_6', 'Gene_10', 'Gene_12', 'Gene_20'])


    if not method in ["leiden", "louvain"]:
        raise ValueError("method argument should be leiden or louvain")
    if type(celltype) == str:
        cluster_subset = _subset_adata(adata, adata.obs.get(celltype_label) == celltype)
    elif type(celltype) == tuple:
        print('hey I am a tuple', celltype_label)
        filter = adata.obs.get(celltype_label) == "NONE"
        print(filter)
        for i in range(len(celltype)):
            filter = filter | (adata.obs.get(celltype_label) == celltype[i])
        cluster_subset = _subset_adata(adata, filter)
        print(cluster_subset.var, filter)
    else:
        sys.exit("specify cluster input as a string or tuple")

    cluster_subset.raw = cluster_subset

    # identify highly variable genes
    sc_highly_variable_genes(
        cluster_subset,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        inplace=True,
        batch_key=batch_key,
    )

    if batch_key != None:
        hvglist = cluster_subset.var["highly_variable"].copy()
        hvglist.loc[
            cluster_subset.var["highly_variable_nbatches"]
            >= len(set(cluster_subset.obs[batch_key])) / n_shared,
        ] = True
        cluster_subset.var["highly_variable"] = hvglist.copy()

    if show_plot_filter:
        pl_highly_variable_genes(cluster_subset, show=True)
    logging.info(
        "In total",
        str(sum(cluster_subset.var.highly_variable)),
        "highly variable genes selected within cluster",
    )

    # apply filter
    cluster_subset = _subset_adata(
        cluster_subset, cluster_subset.var.highly_variable, axis=1, raw=False
    )

    # perform further processing
    # log1p(cluster_subset) # data already logged
    if regress_out_key is not None:
        regress_out(cluster_subset, keys=regress_out_key)
    sc_scale(cluster_subset, max_value=10)
    sc_pca(
        cluster_subset, random_state=random_seed, svd_solver="arpack"
    )  # using `svd_solver='arpack' ensures that the PCA leads to reproducible results
    neighbors(cluster_subset, n_neighbors=10, random_state=random_seed)
    umap(cluster_subset, random_state=random_seed)
    if method == "louvain":
        louvain(cluster_subset, resolution=resolution, random_state=random_seed)
    if method == "leiden":
        leiden(cluster_subset, resolution=resolution, random_state=random_seed)

    return cluster_subset """
""" 
adata = simulated_pbmc3k_processed()
#print(adata)
adata_subset = recluster(adata, celltype=('0', '1', '3', '6'), resolution = 1.3)
#print(adata_subset)

#print('-----------------------------------------')
#print('-----------------------------------------')
#print('-----------------------------------------')

base_adata = bc.datasets.pbmc3k_processed()

base_adata_subset = recluster(base_adata, celltype=('0', '1', '3', '6'), resolution = 1.3)
#print(base_adata_subset)
#print(base_adata_subset.var)
#t = sc.pl.umap(base_adata_subset, color = ['leiden', 'CD3G', 'CD8A', 'CD4', 'IL7R', 'NKG7', 'GNLY'])
#print(t)
#print(adata_subset.var, adata.var)
sc.pl.umap(adata_subset, color = ['leiden', 'Gene_4', 'Gene_5', 'Gene_6', 'Gene_10', 'Gene_12', 'Gene_20'])
"python -c 'from _helper import *; print(simulated_pbmc3k_raw())'" """