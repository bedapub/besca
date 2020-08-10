"""Builtin Datasets.
"""

import os
import urllib.request
from urllib.error import URLError

import pkg_resources

import pandas
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
        adata = read(filename, backup_url = url, cache= False)
    except Exception:
        raise URLError(f'\n\n\n {filename} could not be downloaded from {url}; \n Please download it manually and store it in your besca installation: besca/datasets/data/')
    return adata




def Baron2016_raw():
    """Raw counts from GSE84133 Pancreatic Islets Dataset [Baron M, Veres A, Wolock SL, et al. A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure. Cell Syst. 2016]

    The data consists of raw gene expression counts of single cells from pancreatic islets.

    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Baron2016_raw()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/baron_raw.h5ad')
    adata = check_dl( filename, url ='https://zenodo.org/record/3968315/files/baron2016_raw.h5ad?download=1')
    return adata


def Baron2016_processed():
    """Processed data from GSE84133 Pancreatic Islets Dataset [Baron M, Veres A, Wolock SL, et al. A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure. Cell Syst. 2016]

    The data consists of processed single cell expression data from pancreatic islets. Data was filtered and celltypes were annotated.

    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Baron2016_processed()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/baron_processed.h5ad')
    adata = check_dl( filename, url = 'https://zenodo.org/record/3968315/files/baron2016_processed.h5ad?download=1')
    return adata


def Granja2019_citeSeq():
    """Citeseq table from GSE139369 mixed-phenotype acute leukemia; healthy sample only
      [Granja JM, Klemm S, McGinnis LM, Kathiria AS et al. Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype acute leukemia. Nat Biotechnol 2019 Dec;37(12):1458-1465. PMID: 31792411.]
    
    Returns
    -------
    table : panda dataframe 
            
    Example
    -------

    >>> import besca as bc
    >>> protein_table = bc.datasets.Granja2019_citeSeq()
    >>> protein_table
    
    """    
    filename = pkg_resources.resource_filename('besca', 'datasets/data/scADT_Healthy_counts.csv')
    urllib.request.urlretrieve ("https://zenodo.org/record/3944753/files/scADT_Healthy_counts.csv?download=1", filename)
    citeSeq  = pandas.read_csv(filename, sep = '\t')
    return citeSeq





def Granja2019_raw():
    """Raw data from Granja et al. Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype acute leukemia. Nat Biotechnol. 2019.
    
    The data consists of raw gene expression and CITE-seq counts of single cells from healthy bone marrow, CD34+ bone marow, peripheral blood, and MPAL donors.
    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Raw data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Granja2019_raw()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/Granja2019_raw.h5ad')
    adata = check_dl( filename, url = 'https://zenodo.org/record/3944753/files/Granja2019_raw.h5ad?download=1')
    return adata



def Granja2019_processed():
    """Processed data from Granja et al. Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype acute leukemia. Nat Biotechnol. 2019.
    
     The data consists of processed gene expression and CITE-seq counts of single cells from healthy bone marrow, CD34+ bone marow, peripheral blood, and MPAL donors. Data was filtered and celltypes were annotated.
    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Granja2019_processed()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/Granja2019.annotated.h5ad')
    adata = check_dl( filename, url = 'https://zenodo.org/record/3944753/files/Granja2019.annotated.h5ad?download=1')
    return adata




def Haber2017_raw():
    """Raw counts from GSE92332 mouse small intestine single cell transcriptomics dataset [Haber AL, Biton M, Rogel N, et al. A single-cell survey of the small intestinal epithelium. Nature. 2017]

    The data consists of raw gene expression counts of single cells from mouse small intestine.

    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Haber2017_raw()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/haber_raw.h5ad')
    adata = check_dl( filename, url ='https://zenodo.org/record/3935782/files/haber_raw.h5ad?download=1')
    return adata


def Haber2017_processed():
    """Processed data from GSE92332 mouse small intestine single cell transcriptomics dataset [Haber AL, Biton M, Rogel N, et al. A single-cell survey of the small intestinal epithelium. Nature. 2017]

    The data consists of processed single cell expression data from mouse small intestine. Data was filtered and celltypes were annotated.

    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Haber2017_processed()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/haber_processed.h5ad')
    adata = check_dl( filename, url = 'https://zenodo.org/record/3935782/files/haber_processed.h5ad?download=1')
    return adata



def Kotliarov2020_raw():
    """Raw counts from Kotliarov et al. Broad immune activation underlies shared set point signatures for vaccine responsiveness in healthy individuals and disease activity in patients with lupus. Nat Med. 2020

    The data consists of raw gene expression and CITE-seq counts of single cells from 53,201 single cells from healthy high and low influenza vaccination responders.
    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Raw data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Kotliarov2020_raw()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/Kotliarov2020_raw.h5ad')
    adata = check_dl( filename, url ='https://zenodo.org/record/3938290/files/Kotliarov2020_raw.h5ad?download=1')
    return adata



def Kotliarov2020_processed():
    """Processed data from Kotliarov et al. Broad immune activation underlies shared set point signatures for vaccine responsiveness in healthy individuals and disease activity in patients with lupus. Nat Med. 2020

    The data consists of processed single cell expression and CITE-seq data from healthy high and low influenza vaccination responders. Data was filtered and celltypes were annotated.
    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Kotliarov2020_processed()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/Kotliarov2020_processed_citeseq_merged_annotated.h5ad')
    adata = check_dl( filename, url = 'https://zenodo.org/record/3938290/files/Kotliarov2020_processed_citeseq_merged_annotated.h5ad?download=1')
    return adata



def Lee2020_raw():
    """Raw counts from GSE132465 Primary Colorectal Cancer Dataset [Lee HO, Hong Y, Etlioglu HE, et al. Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer. Nat Genet. 2020]

    The data consists of raw gene expression counts of single cells from primary colorectal cancer.

    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Lee2020_raw()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/lee2020_raw.h5ad')
    adata = check_dl( filename, url ='https://zenodo.org/record/3967538/files/lee2020_raw.h5ad?download=1')
    return adata


def Lee2020_processed():
    """Processed data from GSE132465 Primary Colorectal Cancer Dataset [Lee HO, Hong Y, Etlioglu HE, et al. Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer. Nat Genet. 2020]

    The data consists of processed single cell expression data from primary colorectal cancer. Data was filtered and celltypes were annotated.

    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Lee2020_processed()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/lee2020_processed.h5ad')
    adata = check_dl( filename, url = 'https://zenodo.org/record/3967538/files/lee2020_processed.h5ad?download=1')
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
    adata =  check_dl( filename, url = 'https://zenodo.org/record/3948150/files/pbmc3k_raw.h5ad?download=1')
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
    adata = check_dl( filename, url='https://zenodo.org/record/3948150/files/pbmc3k_filtered.h5ad?download=1')
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
    adata = check_dl(filename, url='https://zenodo.org/record/3948150/files/pbmc3k_processed.h5ad?download=1')
    return adata



def Peng2019_raw():
    """Raw counts from Peng2019 PRJCA001063 Pancreatic Ductal Adenocarcinoma Dataset [Peng J, Sun BF, Chen CY, et al. Single-cell RNA-seq highlights intra-tumoral heterogeneity and malignant progression in pancreatic ductal adenocarcinoma [published correction appears in Cell Res. 2019]

    The data consists of raw gene expression counts of single cells from Pancreatic Ductal Adenocarcinoma.

    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Peng2019_raw()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/StdWf1_PRJCA001063_CRC_besca2.raw.h5ad')
    adata = check_dl( filename, url ='https://zenodo.org/record/3969339/files/StdWf1_PRJCA001063_CRC_besca2.raw.h5ad?download=1')
    return adata


def Peng2019_processed():
    """Processed data from Peng2019 PRJCA001063 Pancreatic Ductal Adenocarcinoma Dataset [Peng J, Sun BF, Chen CY, et al. Single-cell RNA-seq highlights intra-tumoral heterogeneity and malignant progression in pancreatic ductal adenocarcinoma [published correction appears in Cell Res. 2019]

    The data consists of processed single cell expression data from Pancreatic Ductal Adenocarcinoma. Data was filtered and celltypes were annotated.
    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Peng2019_processed()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad')
    adata = check_dl( filename, url = 'https://zenodo.org/record/3969339/files/StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad?download=1')
    return adata



def Segerstolpe2016_processed():
    """Processed data from Segerstolpe et al. Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and Type 2 Diabetes. Cell Metab. 2016.
    
    The data consists of processed gene expression from thousands of human islet cells from healthy and type 2 diabetic human donors. Data was filtered and celltypes were annotated.
    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
            
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.Segerstolpe2016_processed()
    >>> adata
    
    """

    filename = pkg_resources.resource_filename('besca', 'datasets/data/Segerstolpe2016_processed.h5ad')
    adata = check_dl(filename, url = 'https://zenodo.org/record/3928276/files/Segerstolpe2016_processed.h5ad?download=1')
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
    adata = check_dl(filename, url='https://zenodo.org/record/3960617/files/Smillie2019_raw.h5ad?download=1')
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
    adata = check_dl(filename, url='https://zenodo.org/record/3960617/files/Smillie2019_processed.h5ad?download=1')
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
