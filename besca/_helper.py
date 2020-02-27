from anndata import AnnData
import mygene
import sys
from pandas import DataFrame

def subset_adata(adata,
                 filter_criteria,
                 raw = True,
                 axis = 0):
    """ Subset AnnData object into new object

    parameters
    ----------

    adata: `AnnData`
        complete AnnData object

    raw: `bool` | default = True
        boolian indicator if the subset should be initialized with the data contained in adata.raw or not

    axis: `0` or `1` | default = 0
        if 0 then then the filter criteria is applied as a selector for rows, if 1 then the filter criteria is applied as a selector for columns 

    returns
    -------

    AnnData
        if raw = True the AnnData subset initialized with the adata.raw otherwise the AnnData subset initialized with adata

    """
    if axis == 0:
        if len(filter_criteria) == adata.shape[0]:
            # RMK : SLICING HERE WITH ":"" lead to a bug, which should be corected in next anndata version (> 0.22)).
            # Please check before correcting this formulation.
            subset_full = adata[filter_criteria, ].copy()
        else:
            sys.exit('filter criteria used has a different length than the number of rows in AnnData object')
    if axis == 1:
        if len(filter_criteria) == adata.shape[1]:
            subset_full = adata[:, filter_criteria].copy()
        else:
            sys.exit('filter criteria used has a different length than the number of columns in AnnData object')    
    if raw:
        subset = AnnData(X = subset_full.raw.X)
        subset.var_names = subset_full.raw.var_names
        subset.obs_names = subset_full.obs_names
        subset.obs = subset_full.obs
        subset.var = subset_full.raw.var

        return(subset)

    else:
        return(subset_full)

def convert_ensembl_to_symbol(gene_list,
                              species = 'human'):
    """Convert ENSEMBL gene ids to SYMBOLS
    Uses the python package mygene to look up the supplied list of ENSEMBLE Ids and return
    the equivalent list of Symbols. Species needs to be supplied.

    parameters
    ----------
    gene_list: `list`
        list of ensemble_Ids that need to be converted
    species: `str` | default = 'human'
        string identifying the species of the supplied Ensemble IDs

    returns
    -------
    list
        List containing the converted Symbols

    """

    mg = mygene.MyGeneInfo()
    gene_symbols = mg.querymany(gene_list , scopes='ensembl.gene', fields='symbol', species=species)

    symbols = []
    for x in gene_symbols:
        symbols.append(x.get('symbol'))

    return(symbols)

def convert_symbol_to_ensembl(gene_list,
                              species = 'human'):
    """Convert SYMBOLS to ENSEMBL gene ids
    Uses the python package mygene to look up the supplied list of SYMBOLS and return
    the equivalent list of ENSEMBLE GENEIDs. Species needs to be supplied.

    Note: this can result in an error when non-unqiue Symbols are supplied.

    parameters
    ----------
    gene_list: `list`
        list of ensemble_Ids that need to be converted
    species: `str` | default = 'human'
        string identifying the species of the supplied Ensemble IDs

    returns
    -------
    list
        List containing the converted ENSEMBLE gene ids

    """

    mg = mygene.MyGeneInfo()
    gene_symbols = mg.querymany(gene_list , scopes='symbol', fields='ensembl.gene', species=species)

    ensembl = []
    for x in gene_symbols:
        ensembl.append(x.get('ensembl'))

    return(ensembl)

def get_raw(adata):
    """Extract the AnnData object saved in adata.raw

    Returns the AnnData object saved in adata.raw as a normal AnnData object with which you can continue working.
    This includes all of the additional annotation.

    parameters
    ----------
    adata: `AnnData`
        complete AnnData object

    returns
    -------
    AnnData
        returns the AnnData object contained in .raw with all relevant annotation

    Examples
    --------

    >>> import besca as bc
    >>> adata = bc.datasets.pbmc3k_processed()
    >>> print(adata)
    >>> adata_raw = bc.get_raw(adata)
    >>> print(adata_raw)

    """

    if adata.raw is None:
        sys.exit('please pass an AnnData object that contains .raw')

    #initialize new AnnData object
    adata_raw = AnnData(X = adata.raw.X)

    #add varibale attributes from adata.raw
    adata_raw.var_names = adata.raw.var_names
    adata_raw.var = adata.raw.var

    #add observations (same as in adata!)
    adata_raw.obs_names = adata.obs_names
    adata_raw.obs = adata.obs

    #return new adata object
    return(adata_raw)


def get_means(adata,mycat):
    """ Calculates average and fraction expression per category in adata.obs
    Based on an AnnData object and an annotation category (e.g. louvain) returns 
    average expression and fraction cells expressing gene per category
    parameters
    ----------
    adata: AnnData
      an AnnData object
    mycat: str
      the category for summarisation (e.g. louvain, cell_names)
    returns
    -------
    average_obs
        average gene expression per category
    fraction_obs
        fraction cells expressing a gene per category
    """
    gene_ids = adata.raw.var.index.values
    try:
        x = adata.obs[mycat]
        adata.obs[mycat]=adata.obs[mycat].astype('category')
        clusters = adata.obs[mycat].cat.categories
        obs = adata.raw[:,gene_ids].X.toarray()
        obs = DataFrame(obs,columns=gene_ids,index=adata.obs[mycat])
        average_obs = obs.groupby(level=0).mean()
        obs_bool = obs.astype(bool)
        fraction_obs = obs_bool.groupby(level=0).sum()/obs_bool.groupby(level=0).count()
    except KeyError:
        print("Oops!  The adata object does not have the specified column. Options are: ")
        print (list(adata.obs.columns))
        average_obs=None
        fraction_obs=None
    return(average_obs, fraction_obs)
