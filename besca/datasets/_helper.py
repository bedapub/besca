"""Simulated datasets.
"""

from anndata import AnnData
import numpy as np
from scipy.sparse import csr_matrix
import random
from pandas.api.types import CategoricalDtype


def simulated_Kotliarov2020_processed(n_var = 25, n_obs = 250):
    """Simulated Kotliarov2020 processed data (used for doctests)

    The real Kotliarov2020 processed data can be found under besca.datasets.Kotliarov2020_processed()

    Parameters
    ----------
    n_var: `int`
        optional, n_var size of AnnData matrix, default value 25
    n_obs: `int`
        optional, n_obs size of AnnData matrix, default value 250
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    """
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


def simulated_pbmc3k_raw(n_var = 100, n_obs = 10):
    """Simulated 3k PBMC raw data (used for doctests)

    The real 3k PBMC raw data can be found under besca.datasets.pbmc3k_raw()

    Parameters
    ----------
    n_var: `int`
        optional, n_var size of AnnData matrix, default value 100
    n_obs: `int`
        optional, n_obs size of AnnData matrix, default value 10
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    """
    counts = csr_matrix(np.random.poisson(1, size=(n_obs, n_var)), dtype=np.float32)
    adata = AnnData(counts)
    obs_names = [''.join((random.choice('ACGT') for x in range(14))) for i in range(adata.n_obs)]
    obs_names = [''.join(f"{name}-1") for name in obs_names]
    adata.obs_names = obs_names
    adata.obs["CELL"] = obs_names
    adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
    adata.var["SYMBOL"] = [f"Gene_{i:d}" for i in range(adata.n_vars)]
    return adata


def simulated_pbmc3k_processed(n_var = 25, n_obs = 250):
    """Simulated 3k PBMC processed data (used for doctests)

    The real 3k PBMC processed data can be found under besca.datasets.pbmc3k_processed()

    Parameters
    ----------
    n_var: `int`
        optional, n_var size of AnnData matrix, default value 25
    n_obs: `int`
        optional, n_obs size of AnnData matrix, default value 250
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    """
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
    adata.obs["dblabel"] = adata.obs["dblabel"].astype("category")
    adata.obs["celltype3"] = np.random.choice(["naive B cell", "CD1c-positive myeloid dendritic cell", "non-classical monocyte"], size=(adata.n_obs,))
    adata.obs["celltype3"] = adata.obs["celltype3"].astype("category")
    adata.obs["leiden"] =  [str(random.randint(0,12)) for x in range(adata.n_obs)]
    status_type = CategoricalDtype(categories=[f"{x}" for x in range(13)], ordered=True)
    adata.obs["leiden"] = adata.obs["leiden"].astype(status_type)
    adata.obsm["X_umap"] = np.random.normal(0, 1, size=(adata.n_obs, 2))
    return adata


def simulated_Haber2017_processed(n_var = 10, n_obs = 1000):
    """Simulated Haber2017 processed data (used for doctests)

    The real Haber2017 processed data can be found under besca.datasets.Haber2017_processed()

    Parameters
    ----------
    n_var: `int`
        optional, n_var size of AnnData matrix, default value 10
    n_obs: `int`
        optional, n_obs size of AnnData matrix, default value 1000
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    """
    raw_data = [[y if bool(random.getrandbits(1)) else -abs(y) for y in x] for x in np.random.exponential(1.6, size=(n_obs, n_var))]
    counts = csr_matrix(raw_data, dtype=np.float32)
    raw_adata = AnnData(counts)
    obs_names = [''.join((random.choice('ACGT') for x in range(14))) for i in range(raw_adata.n_obs)]
    obs_names = [''.join(f"Donor1.{name}") for name in obs_names]
    raw_adata.obs["CELL"] = obs_names
    raw_adata.var_names = [f"Gene_{i:d}" for i in range(raw_adata.n_vars)]
    raw_adata.var["SYMBOL"] = [f"Gene_{i:d}" for i in range(raw_adata.n_vars)]
    adata = AnnData(np.asarray(raw_data))
    adata.raw = raw_adata
    adata.obs["CELL"] = obs_names
    adata.obs["leiden"] =  [str(random.randint(0,38)) for x in range(adata.n_obs)]
    status_type = CategoricalDtype(categories=[f"{x}" for x in range(39)], ordered=True)
    adata.obs["leiden"] = adata.obs["leiden"].astype(status_type)
    adata.obs["donor"] =  [''.join(f"donor{random.randint(1,2)}") for x in range(adata.n_obs)]
    adata.obs["donor"] = adata.obs["donor"].astype("category")
    adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
    adata.var["SYMBOL"] = [f"Gene_{i:d}" for i in range(adata.n_vars)]
    return adata


def simulated_Baron2016_processed(n_var = 100, n_obs = 20):
    """Simulated Baron2016 processed data (used for doctests)

    The real Baron2016 processed data can be found under besca.datasets.Baron2016_processed()

    Parameters
    ----------
    n_var: `int`
        optional, n_var size of AnnData matrix, default value 100
    n_obs: `int`
        optional, n_obs size of AnnData matrix, default value 20
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    """
    raw_data = np.random.exponential(1.6, size=(n_obs, n_var))
    counts = csr_matrix(raw_data, dtype=np.float32)
    raw_adata = AnnData(counts)
    raw_adata.var_names = [f"Gene_{i:d}" for i in range(raw_adata.n_vars)]
    raw_adata.var["SYMBOL"] = [f"Gene_{i:d}" for i in range(raw_adata.n_vars)]
    adata = AnnData(np.asarray(raw_data))
    adata.raw = raw_adata
    adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
    adata.var["SYMBOL"] = [f"Gene_{i:d}" for i in range(adata.n_vars)]
    adata.obs["celltype2"] = np.random.choice(["pancreatic acinar cell", "pancreatic D cell", "mast cell", "myeloid leukocyte"], size=(adata.n_obs,))
    adata.obs["celltype2"] = adata.obs["celltype2"].astype("category")
    adata.obs["assigned_cluster"] = np.random.choice(["acinar", "beta", "alpha", "t_cell", "quiescent_stellate"], size=(adata.n_obs,))
    adata.obs["assigned_cluster"] = adata.obs["assigned_cluster"].astype("category")
    return adata
