#this contains wrapper functions to export data into the standard format for the standard pipeline

#import other functions
from ..export._export import X_to_mtx, raw_to_mtx, louvain, clustering,  labeling, analysis_metadata, ranked_genes, labeling_info
from ..import._read import read_mtx
from .. import _logging as logs
from os.path import join
from os import makedirs
import os
from time import time
import logging
from ..pp._filtering import filter

def export_cp10k(adata, basepath):
    """ Export raw cp10k to FAIR format for loading into database

    wrapper function for X_to_mtx with correct folder structure for loading into database.

    parameters
    ----------
    adata: `ÀnnData`
        AnnData object that is to be exported
    basepath: `str`
        root path to the Analysis folder (i.e. ../analyzed/<ANALYSIS_NAME>)

    returns
    -------
    None
        writes to file

    """
    
    #call wrapper function
    X_to_mtx(adata=adata,
             outpath=join(basepath, 'normalized_counts', 'cp10k'),
             write_metadata=True,
             geneannotation='SYMBOL',
             additional_geneannotation='ENSEMBL')


def export_regressedOut(adata, basepath):
    """ Export regressedOut to FAIR format for loading into database

    wrapper function for X_to_mtx with correct folder structure for loading into database.

    parameters
    ----------
    adata: `ÀnnData`
        AnnData object that is to be exported
    basepath: `str`
        root path to the Analysis folder (i.e. ../analyzed/<ANALYSIS_NAME>)

    returns
    -------
    None
        writes to file
        
    """
    X_to_mtx(adata, 
             outpath=join(basepath, 'normalized_counts', 'regressedOut'), 
             geneannotation='SYMBOL',
             write_metadata= True,
             additional_geneannotation='ENSEMBL')



## LEFT FOR COMPATIBILTIY SEE export_clustering for new function
def export_clustering(adata, basepath, method):
    """ Export louvain to cell mapping to FAIR format for loading into database

    wrapper function for louvain and labeling_info with correct folder structure/names
    for loading into the database.

    parameters
    ----------
    adata: `ÀnnData`
        AnnData object that is to be exported
    basepath: `str`
        root path to the Analysis folder (i.e. ../analyzed/<ANALYSIS_NAME>)
    method: `str`
        method of clustering used previously, should be leiden or louvain

    returns
    -------
    None
        writes to file
    """  

    clustering(adata, 
         outpath = join(basepath, 'labelings', method), method = method)
    labeling_info(outpath=join(basepath, 'labelings', method), description=method+' clustering', method=method)



def export_metadata(adata,
                    basepath,
                    n_pcs = 3,
                    umap = True,
                    tsne = False):
    """ Export metadata in FAIR format for loading into database

    wrapper function for analysis_metadata with correct folder structure/names
    for loading into the database.

    parameters
    ----------
    adata: `ÀnnData`
        AnnData object that is to be exported
    basepath: `str`
        root path to the Analysis folder (i.e. ../analyzed/<ANALYSIS_NAME>)
    n_pcs: `int` | default = 3
        number of PCA components to export
    umap: `bool` | default = True
        boolian indicator if umap coordinates should be exported
    tsne: `bool` | default = False
        boolian indicator if tSNE coordinates should be exported

    returns
    -------
    None
        writes to file
    
    """   
    
    analysis_metadata(adata, 
                      outpath=join(basepath), 
                      n_pcs= n_pcs, 
                      umap=umap, 
                      tsne=tsne)


def export_rank(adata, basepath, type = 'wilcox', labeling_name = 'louvain'):
    """ Export ranked genes to FAIR format for loading into database

    wrapper function for ranked_genes with correct folder structure/names
    for loading into the database.

    parameters
    ----------
    adata: `ÀnnData`
        AnnData object that is to be exported
    basepath: `str`
        root path to the Analysis folder (i.e. ../analyzed/<ANALYSIS_NAME>)
    type: `str` | default = 'wilcox'
        indicator of the statistical method employed, can be one of: 'wilcox' or 't-test overest var'  or 't-test'
    labelingname: `str` | default = `louvain`
        labeling that will be exported

    returns
    -------
    None
        writes to file
    
    """   
    #export rank files
    ranked_genes(adata=adata, 
                 outpath=join(basepath, 'labelings', labeling_name),
                 type= type)

def export_celltype(adata, basepath):
    """ Export celltype annotation to cell mapping in FAIR format for loading into database

    wrapper function for labeling and labeling_info with correct folder structure/names
    for loading into the database.

    parameters
    ----------
    adata: `ÀnnData`
        AnnData object that is to be exported
    basepath: `str`
        root path to the Analysis folder (i.e. ../analyzed/<ANALYSIS_NAME>)

    returns
    -------
    None
        writes to file
    
    """   

    #export celltypes
    labeling(adata=adata,
             outpath = join(basepath, 'labelings', 'celltype'), 
             column='celltype')

    labeling_info(outpath = join(basepath, 'labelings', 'celltype'),
                  method='manual celltype annotation based on marker expression',
                  annotated_version_of='louvain', 
                  expert=True, 
                  default=False, 
                  public=False, 
                  reference=True, 
                  description='manual celltype annotation based on the expression of marker genes')
