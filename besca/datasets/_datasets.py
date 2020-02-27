"""Builtin Datasets.
"""

import os
from scanpy.api import read

def pbmc3k_raw():
    """3k PBMCs from 10x Genomics raw.

    The data consists of 3k PBMCs from a Healthy Donor and is freely available
    from 10x Genomics (`here
    <http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz>`__
    from this `webpage
    <https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k>`__).
    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.

    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.pbmc3k_raw()
    >>> adata

    """

    filename = os.path.dirname(__file__) + '/data/pbmc3k_raw.h5ad'
    adata = read(filename, cache = True)
    return adata

def pbmc3k_filtered():
    """3k PBMCs from 10x Genomics filtered.

    The data consists of filtered data of 3k PBMCs from a Healthy Donor.
    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.

    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.pbmc3k_filtered()
    >>> adata

    """

    filename = os.path.dirname(__file__) + '/data/pbmc3k_filtered.h5ad'
    adata = read(filename, cache = True)
    return adata

def pbmc3k_processed():
    """3k PBMCs from 10x Genomics processed

    The data consists of filtered data of 3k PBMCs from a Healthy Donor that has been
    completely processed (i.e. PCA, nearest neighbors, UMAP, tSNE and rank_genes_groups().
    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.pbmc3k_processed()
    >>> adata
    
    """

    filename = os.path.dirname(__file__) + '/data/pbmc3k_processed.h5ad'
    adata = read(filename, cache = True)
    return adata

def pbmc_storage_raw():
    """PBMCs at 3 storage conditions raw

    The data consists of raw data of PBMCs from 3 healthy donor that were stored under 3 different conditions: fresh, frozen, 24_RT.
    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.pbmc_storage_raw()
    >>> adata
    
    """

    filename = os.path.dirname(__file__) + '/data/pbmc_storage_raw_downsampled.h5ad'
    adata = read(filename, cache=True)
    return adata
  

def pbmc_storage_processed():
    """PBMCs at 3 storage conditions raw

    The data consists of raw data of PBMCs from 3 healthy donor that were stored under 3 different conditions: fresh, frozen, 24_RT.
    Data was filtered, batch corrected using BBKNN and celltypes were annotated.
    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.

    Example
    -------
    >>> import besca as bc
    >>> adata = bc.datasets.pbmc_storage_processed()
    >>> adata
    """
    filename = os.path.dirname(__file__) + '/data/pbmc_storage_processed_downsampled.h5ad'
    adata = read(filename, cache=True)
    return adata



def load_immune_signatures(refined=True):
    """ This function load immune signatures provided 
    within BESCA.

    Parameters
    ----------
    refined: `logical` | default = True
        if True, a large set of gene signatures including subpopulations 
        will be provided.
        if False only large markers will be returned.

    Returns
    -------
    dict : a dictionnary containing a key the name of the 
    signature; as value an array of gene symbols
            
    Example
    -------

    >>> import besca as bc
    >>> immune_sig = bc.datasets.load_immune_signature()
    >>> immune_sig
    """
    if refined:
        filename = os.path.dirname(__file__) + '/genesets/HumanCD45p_scseqCMs6.gmt'
    else: 
        filename = os.path.dirname(__file__) + '/genesets/Immune.gmt'
    file = open(filename, 'r')
    Lines = file.readlines()
    mymarkers = {}
    for line in Lines:
        ll = line.strip().split('\t')
        mymarkers[ll[0]] = ll[2:len(ll)]
    file.close()
    return mymarkers
