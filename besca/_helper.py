from anndata import AnnData
import mygene
import sys
from pandas import DataFrame
from scipy import sparse
from pandas import concat


def subset_adata(adata, filter_criteria, raw=True, axis=0):
    """Subset AnnData object into new object

    parameters
    ----------

    adata: `AnnData`
        complete AnnData object
    filter_criteria: bool serie
        filtering requirement along the axis chosen.
    raw: `bool` | default = True
        boolian indicator if the subset should be initialized with the data contained in adata.raw or not

    axis: `0` or `1` | default = 0
        if 0 then then the filter criteria is applied as a selector for rows, if 1 then the filter criteria is applied as a selector for columns

    returns
    -------

    AnnData
        if raw = True the AnnData subset initialized with the adata.raw otherwise the AnnData subset initialized with adata


    Examples
    --------

    >>> import besca as bc
    >>> adata = bc.datasets.pbmc3k_processed()
    >>> adata_BCELL =  bc.subset_adata( adata, filter_criteria= adata.obs.get('dblabel') == 'naive B cell')

    """
    if axis == 0:
        if len(filter_criteria) == adata.shape[0]:
            # RMK : SLICING HERE WITH ":"" lead to a bug, which should be corected in next anndata version (> 0.22)).
            # Please check before correcting this formulation.
            subset_full = adata[
                filter_criteria,
            ].copy()
        else:
            sys.exit(
                "filter criteria used has a different length than the number of rows in AnnData object"
            )
    if axis == 1:
        if len(filter_criteria) == adata.shape[1]:
            subset_full = adata[:, filter_criteria].copy()
        else:
            sys.exit(
                "filter criteria used has a different length than the number of columns in AnnData object"
            )
    if raw:
        subset = AnnData(X=subset_full.raw.X)
        subset.var_names = subset_full.raw.var_names
        subset.obs_names = subset_full.obs_names
        subset.obs = subset_full.obs
        subset.var = subset_full.raw.var

        return subset

    else:
        return subset_full


def convert_ensembl_to_symbol(gene_list, species="human"):
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
    gene_symbols = mg.querymany(
        gene_list, scopes="ensembl.gene", fields="symbol", species=species
    )

    symbols = []
    for x in gene_symbols:
        symbols.append(x.get("symbol"))

    return symbols


def convert_symbol_to_ensembl(gene_list, species="human"):
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
    gene_symbols = mg.querymany(
        gene_list, scopes="symbol", fields="ensembl.gene", species=species
    )

    ensembl = []
    for x in gene_symbols:
        ensembl.append(x.get("ensembl"))

    return ensembl


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
        sys.exit("please pass an AnnData object that contains .raw")

    # initialize new AnnData object
    adata_raw = AnnData(X=adata.raw.X)

    # add varibale attributes from adata.raw
    adata_raw.var_names = adata.raw.var_names
    adata_raw.var = adata.raw.var

    # add observations (same as in adata!)
    adata_raw.obs_names = adata.obs_names
    adata_raw.obs = adata.obs

    # return new adata object
    return adata_raw


def get_means(adata, mycat):
    """Calculates average and fraction expression per category in adata.obs
    Based artihmetic mean expression and fraction cells expressing gene per category
    (works on linear scale)
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
        adata.obs[mycat] = adata.obs[mycat].astype("category")
        clusters = adata.obs[mycat].cat.categories
        obs = expm1.adata.raw[:, gene_ids].X.toarray()
        obs = DataFrame(obs, columns=gene_ids, index=adata.obs[mycat])
        average_obs = log1p(obs.groupby(level=0).mean())
        obs_bool = obs.astype(bool)
        fraction_obs = (
            obs_bool.groupby(level=0).sum() / obs_bool.groupby(level=0).count()
        )
    except KeyError:
        print(
            "Oops!  The adata object does not have the specified column. Options are: "
        )
        print(list(adata.obs.columns))
        average_obs = None
        fraction_obs = None
    return (average_obs, fraction_obs)


def get_gmeans(adata, mycat):
    """Calculates average and fraction expression per category in adata.obs
    Based on an AnnData object and an annotation category (e.g. louvain) returns
    geometric mean expression and fraction cells expressing gene per category
    (works on log scale)
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
        adata.obs[mycat] = adata.obs[mycat].astype("category")
        clusters = adata.obs[mycat].cat.categories
        obs = adata.raw[:, gene_ids].X.toarray()
        obs = DataFrame(obs, columns=gene_ids, index=adata.obs[mycat])
        average_obs = obs.groupby(level=0).mean()
        obs_bool = obs.astype(bool)
        fraction_obs = (
            obs_bool.groupby(level=0).sum() / obs_bool.groupby(level=0).count()
        )
    except KeyError:
        print(
            "Oops!  The adata object does not have the specified column. Options are: "
        )
        print(list(adata.obs.columns))
        average_obs = None
        fraction_obs = None
    return (average_obs, fraction_obs)


def concate_adata(adata1, adata2):
    """concatenate two adata objects based on the observations

    this function also merges the objects saved in .raw and generates a new combined.raw.
    The obs from adata1 are preserved. Those from adata2 are lost.

    parameters
    ----------
    adata1: ´AnnData´
        first Anndata object that is to be concatenated (the obs and obs_names will be taken from here)
    adata2: ´AnnData´
        second Anndata objec that is to be concatenated

    returns
    -------
    AnnData
        returns the AnnData object contained in .raw with all relevant annotation

    """

    # convert both to sparse matrices
    adata1_X = sparse.csr_matrix(adata1.X)
    adata2_X = sparse.csr_matrix(adata2.X)

    # hstack (i.e. horizontally stack them)
    X = sparse.hstack([adata1_X, adata2_X])

    # convert back to sparse matrix (for some reason hstack returns coo matrix)
    X = sparse.csr_matrix(X)

    # observations remain unchanged so will use obs from adata1
    if sum(adata1.obs.CELL != adata2.obs.CELL) != 0:
        sys.exit(
            "The observations in adata1 and adata2 are not the same. These two anndata objects can not be concatenated."
        )

    obs = adata1.obs
    obs_names = adata1.obs_names

    # vars, var_names are concatenated from the two dataframes
    var = concat([adata1.var, adata2.var])
    var_names = adata1.var_names.tolist() + adata2.var_names.tolist()

    if var.index.tolist() != var_names:
        sys.exit("var and var_names index does not match. Something went wrong.")

    # create our complete combined object
    adata_combined = AnnData(X=X, obs=obs, var=var)

    if any([adata1.raw is not None, adata2.raw is not None]):
        if any([adata1.raw is None, adata2.raw is None]):
            sys.exit("Only one of the anndata objects contains a .raw!")

        # get the raw object
        adata1_raw = get_raw(adata1)
        adata2_raw = get_raw(adata2)

        # check that the raw_matrices are the same
        if adata1_raw.X.shape == adata2_raw.X.shape:
            X_raw = adata1_raw.X
        else:
            # need to merge the raw matrices
            # convert both to sparse matrices
            adata1_raw_X = sparse.csr_matrix(adata1_raw.X)
            adata2_raw_X = sparse.csr_matrix(adata2_raw.X)

            # hstack (i.e. horizontally stack them)
            X_raw = sparse.hstack([adata1_raw_X, adata2_raw_X])

            # convert back to sparse matrix (for some reason hstack returns coo matrix)
            X_raw = sparse.csr_matrix(X_raw)

        # get the raw var and obs names

        # check if .raw.var matches
        if adata1_raw.var.shape == adata2_raw.var.shape:
            var_raw = adata1_raw.var
            var_names = adata1_raw.var_names
        else:
            var_raw = concat([adata1_raw.var, adata2_raw.var])
            var_names_raw = (
                adata1_raw.var_names.tolist() + adata2_raw.var_names.tolist()
            )

            if var_raw.index.tolist() != var_names_raw:
                sys.exit(
                    "var and var_names index does not match in .raw. Something went wrong."
                )

        obs_raw = adata1_raw.obs
        obs_names_raw = adata1_raw.obs_names

        # create our adata_combined_raw object (this will be the same for all of the other combined objects as well)
        adata_combined_raw = AnnData(X=X_raw, obs=obs_raw, var=var_raw)

        # add our new raw objec to our combined adata object
        adata_combined.raw = adata_combined_raw

    return adata_combined
