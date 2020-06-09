"""Builtin Datasets.
"""

import os
from urllib.error import URLError
import pkg_resources
from scanpy.api import read


def check_dl( filename, url):
    """try to obtain dataset while checking url 
      ----------
    filename: `str`
        path and filename to load
    url: `str` 
        backup url.
            
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    """
    try:
        adata = read(filename, backup_url = url, cache= True)
    except Exception:
        raise URLError(f'\n\n\n {filename} could not be downloaded from {url}; \n Please download it manually and store it in your besca installation: besca/datasets/data/')
    return adata


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

    filename = pkg_resources.resource_filename('besca', 'datasets/data/pbmc3k_raw.h5ad')
    adata =  check_dl( filename, url = 'https://zenodo.org/record/3886414/files/pbmc3k_raw.h5ad?download=1')
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

    filename = pkg_resources.resource_filename('besca', 'datasets/data/pbmc3k_filtered.h5ad')
    adata = check_dl( filename, url='https://zenodo.org/record/3886414/files/pbmc3k_filtered.h5ad?download=1')
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

    filename = pkg_resources.resource_filename('besca', 'datasets/data/pbmc3k_processed.h5ad')
    adata = check_dl(filename, url='https://zenodo.org/record/3886414/files/pbmc3k_processed.h5ad?download=1')
    return adata




def Smillie2019_raw():
    """Raw counts from Smillie et al. Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis. Cell. 2019

    The data consists of raw gene expression counts of single cells from colon mucosa from 7 ulcerative colitis (UC) patients and 10 healthy controls, paired samples (inlamed, non-inflamed for UC, location-matched for healthy): 34 samples. Epithelial (EPI) and lamina propria (LP) fractions enriched in a two-step digestion process. 
    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Smillie2019_raw()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/Smillie2019_raw.h5ad')
    adata = read(filename, cache=True)
    return adata


def Smillie2019_processed():
    """Processed data from Smillie et al. Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis. Cell. 2019

    The data consists of processed single cell expression data from colon mucosa from 7 ulcerative colitis (UC) patients and 10 healthy controls, paired samples (inlamed, non-inflamed for UC, location-matched for healthy): 34 samples. Epithelial (EPI) and lamina propria (LP) fractions enriched in a two-step digestion process. Data was filtered, batch corrected using BBKNN and celltypes were annotated.


    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Smillie2019_processed()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/Smillie2019_processed.h5ad')
    adata = read(filename, cache=True)
    return adata


def Martin2019_raw():
    """Raw counts from Martin et al. Single-Cell Analysis of Crohn's Disease Lesions Identifies a Pathogenic Cellular Module Associated with Resistance to Anti-TNF Therapy. Cell. 2019

    The data consists of raw gene expression counts of single cells from lamina propria cells from inflamed and non-inflamed ileum lesions (and peripheral blood, but not part of GEO dataset) from 11 Crohn’s disease patients: 22 samples

    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Martin2019_raw()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/Martin2019_raw.h5ad')
    adata = check_dl( filename, url ='https://zenodo.org/record/3862132/files/Martin2019_raw.h5ad?download=1')
    return adata


def Martin2019_processed():
    """Processed data from Martin et al. Single-Cell Analysis of Crohn's Disease Lesions Identifies a Pathogenic Cellular Module Associated with Resistance to Anti-TNF Therapy. Cell. 2019

    The data consists of processed single cell expression data from lamina propria cells from inflamed and non-inflamed ileum lesions (and peripheral blood, but not part of GEO dataset) from 11 Crohn’s disease patients: 22 samples. Data was filtered and celltypes were annotated.

    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Martin2019_processed()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/Martin2019_processed.h5ad')
    adata = check_dl( filename, url = 'https://zenodo.org/record/3862132/files/Martin2019_processed.h5ad?download=1')
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
        filename = pkg_resources.resource_filename('besca', 'datasets/genesets/HumanCD45p_scseqCMs6.gmt')
    else: 
        filename = pkg_resources.resource_filename('besca', 'datasets/genesets/Immune.gmt')
    file = open(filename, 'r')
    Lines = file.readlines()
    mymarkers = {}
    for line in Lines:
        ll = line.strip().split('\t')
        mymarkers[ll[0]] = ll[2:len(ll)]
    file.close()
    return mymarkers
