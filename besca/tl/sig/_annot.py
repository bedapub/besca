# this module contains functions for cell type annotation based on signatures in python using scanpy
import numpy as np
import pandas as pd
from scanpy import AnnData
from scipy.stats import mannwhitneyu

def getset(df: pd.DataFrame, signame_complete: str, threshold) -> set:
    """Handles missing signatures aux function for make_anno
    Based on a dataframe of p-values, a signature name and a cutoff check if sign is present
    parameters
    ----------
    df: panda.DataFrame
      a dataframe of p-values per signature and cell cluster
    signame_complete: str
      signature name
    threshold: numpy.float64
      cutoff used for cluster attributions
    returns
    -------
    set
        Set of clusters passing the threshold value
    """

    if bool(df.index.isin([signame_complete]).any()):
        return set(df.columns[df.loc[signame_complete, :] > threshold])
    else:
        return set()


def gene_fr(f: pd.DataFrame, gene_name: str):
    """Returns the fraction of cells expressing a gene in distinct clusters
    Takes as an input a dataframe with fractions per clusters and a gene name
    Returns the fraction genes per cluster
    parameters
    ----------
    f: panda.DataFrame
      a dataframe of fraction genes per cell as outputted by besca
    gene_name: 'str'
      gene name

    returns
    -------
    panda.Series
        a Series of fractions per cluster and gene
    """
    k = f.loc[
        f["Description"].isin(gene_name),
    ].copy()
    k = k.iloc[:, 0 : len(k.columns)]
    return k.mean(axis=0, skipna=True)


def score_mw(f: pd.DataFrame, mymarkers: dict) -> pd.DataFrame:
    """Score Clusters based on a set of immune signatures to generate a df of pvals
    Takes as an input a dataframe with fractions per clusters and a dictionary of signatures
    Performs a Mann-Whitney test per each column and Signature, returns -10logpValues

    Parameters
    ----------
    f: panda.DataFrame
      a dataframe of -log10pvals per signature and cell cluster
    mymarkers: Dictionary
      a dictionary of signatures

    returns
    -------
    panda.DataFrame
        a dataframe of -10logpValues per cluster and signature
    """
    ids = set(f.iloc[:, 2 : len(f.columns)])
    mypFrame = pd.DataFrame(index=mymarkers.keys(), columns=list(ids))
    for key, value in mymarkers.items():
        for i in list(ids):
            mypFrame.loc[key, i] = -10 * np.log(
                mannwhitneyu(
                    x=f.loc[f["Description"].isin(value), :][i],
                    y=f[i],
                    alternative="greater",
                ).pvalue
            )
    return mypFrame


def add_anno(adata: AnnData, cnames: pd.DataFrame, mycol, clusters="leiden"):
    """Adds annotation generated with make_anno to a AnnData object
    Takes as input the AnnData object to which annotation will be appended and the annotation pd
    Generates a pd.Series that can be added as a adata.obs column with annotation at a level

    Parameters
    ----------
    adata: AnnData
      AnnData object that is to be annotated
    cnames: panda.DataFrame
      a list of cluster names generated with make_anno
    mycol: string
      column to be added
    clusters: string
      initial cluster used for annotation (column in adata.obs)

    returns
    -------
    pd.Series
        the annotated cells
    """

    newlab = adata.obs[clusters]
    # mycol='celltype0'

    notmatched = list(set(list(cnames.index)) - set(newlab))
    if len(notmatched) > 0:
        errm = (
            "Please check that you have indexed the correct clustering. "
            + str(len(notmatched))
            + "clusterIDs not found in "
            + clusters
        )
        return errm
    else:
        for i in list(set(newlab)):
            newlab = newlab.replace(i, cnames.loc[i, mycol]).copy()
        return newlab


def read_annotconfig(configfile: str):
    """Reads the configuration file for signature-based hierarhical annotation.

    Parameters
    ----------
     configfile: 'str'
      path to the configuration file, tsv

    returns
    -------
    sigconfig: panda.DataFrame
      a dataframe with the config details
    levsk: list of str
        String of distinct levels based on configuration file
    """
    # read the config file
    sigconfig = pd.read_csv(configfile, sep="\t", index_col=0)
    # Reorder with the specified order. Place better signatures first, only first match will be kept
    sigconfig = sigconfig.sort_values("Order")
    # Consider up to 7 levels
    nochild = list(set(sigconfig.index) - set(sigconfig["Parent"]))
    levs = []
    levs.append(
        list(
            sigconfig.loc[
                sigconfig["Parent"] == "None",
            ].index
        )
    )
    for i in range(0, 7):
        levs.append(
            list(
                sigconfig.loc[
                    sigconfig["Parent"].isin(levs[i]),
                ].index
            )
        )

    levsk = []
    for lev in levs:
        if len(lev) > 0:
            levsk.append(lev)

    return (sigconfig, levsk)


def export_annotconfig(sigconfig, levsk, resultsFolder, analysisName=""):
    """Export the configuration defined in sigconfig and levsk
    Order might changed compared to the original sig. config file.
    But numerical order within a same category (same parents) should be respected.

    Parameters
    ----------
    sigconfig: panda.DataFrame
      a dataframe with the config details
    levsk: list of str
        String of distinct levels based on configuration file
    resultsFolder: 'str'
        results folder. where the annotation will also be stored
    analysisName: 'str'
        will prefix the saved config file
    """
    it_order = 0  # Iterative order.
    config_t = sigconfig.copy()
    get_order = {}
    for child in levsk:
        for element in child:
            get_order[element] = it_order
            it_order += 1
    config_t["Order"] = config_t.index.map(get_order)
    config_t.to_csv(resultsFolder + analysisName + "_config.tsv", sep="\t")


def make_anno(
    df: pd.DataFrame,
    sigscores: dict,
    sigconfig: pd.DataFrame,
    levsk: list,
    lab: str = "celltype",
    toexclude: list = [],
) -> pd.DataFrame:
    """Annotate cell types
    Based on a dataframe of -log10pvals, a cutoff and a signature set generate cell annotation
    Hierarchical model of Immune cell annotation.
    It expects a specific set of signatures with a common prefix (signame) and specified suffixes indicating
    the cell type.

    Parameters
    ----------
    mypFrame: panda.DataFrame
      a dataframe of -log10pvals per signature and cell cluster
    sigscores: dict
      a dictionary with cluster attribution per signature
    sigconfig: panda.DataFrame
      a dataframe with the configuration information
    levsk: list of str
      String of distinct levels based on configuration file
    lab: 'str'
      cell type category base name
    toexclude: list of str
      String of cell types to ignore in the annotation process

    returns
    -------
    cnames: panda.DataFrame
      a dataframe with cluster to cell type attribution per distinct levels
    """

    # make sure that all levels are present in df, remove toexclude
    levskk = []
    # Iterating throuhg the levels and removing excluded ones
    toinc = set(sigscores.keys()) - set(toexclude)
    for x in levsk:
        tmp = [y for y in x if y in toinc]
        levskk.append(tmp)
    # levskk is the updated level list
    sigscoresk = {x: sigscores[x] for x in sigscores.keys() if not x in toexclude}

    # First part, get cluster identities
    myclust = list(df.columns)
    cnames = pd.DataFrame(index=myclust)
    # DATA INITIALISATION; first level
    colname = lab + "0"
    cnames[colname] = None
    for celltype in levskk[0]:
        for ii in myclust:
            if ii in sigscoresk[celltype]:
                if cnames.loc[ii, colname] is None:
                    cnames.loc[ii, colname] = celltype
    # ASSIGNED ALL FIRST LEVEL
    if len(levskk) > 2:
        for j in range(len(levskk) - 1):
            parent_col = lab + str(j)
            colname = lab + str(j + 1)
            cnames[colname] = None
            for ii in myclust:
                # Find child level and removing excluding term
                sublev = set(
                    sigconfig.loc[
                        sigconfig["Parent"] == cnames.loc[ii, parent_col], :
                    ].index
                ).intersection(toinc)
                # getting the right order using levskk
                sublevk = [x for x in levskk[j + 1] if x in sublev]
                assigned_type = False
                for lev in sublevk:
                    if ii in sigscoresk[lev]:
                        if cnames.loc[ii, colname] is None:
                            cnames.loc[ii, colname] = lev
                            assigned_type = True
                # If no value, put parent value
                if not assigned_type:
                    cnames.loc[ii, colname] = cnames.loc[ii, parent_col]
    return cnames


def match_cluster(
    adata: AnnData, obsquery: str, obsqueryval: str, obsref: str = "leiden", cutoff=0.5
) -> list:
    """Matches categories from adata.obs to each other. For a query category specified in obsquery
    and a value specified in obsqueryval, checks which clusters (or other adata.obs categories, obsref)
    contain >50% (or distinct cutoff, cutoff) of cells of the specified kind.

    Parameters
    ----------
    adata: AnnData
      AnnData object
    obsquery: 'str'
      adata.obs category name used for querying
    obsqueryval: 'str'
      adata.obs category name value, present in obsquery
    obsref: 'str'
      adata.obs category name to be returned
    cutoff: 'numpy.float64'
      fraction of positive cells returned
    returns
    -------
    list
        a list of the cluster IDs that match the query label
    """

    myc = []
    myleiden = list(set(adata[adata.obs[obsquery].isin([obsqueryval])].obs[obsref]))
    for i in myleiden:
        mysub = adata[adata.obs[obsref].isin([i])].copy()
        myc.append(len(mysub[mysub.obs[obsquery] == obsqueryval]) / len(mysub))
    return list(np.array(myleiden)[np.array(myc) > cutoff])


def obtain_new_label(
    nomenclature_file: str,
    cnames: pd.DataFrame,
    reference_label: str = "short_dblabel",
    new_label: str = "dblabel",
    new_level: int = None,
) -> pd.DataFrame:

    """Matches the cnames obtained by the make_annot function or a list of label names
    to the db label (standardized label from a nomenclature file).

    parameters
    ---------
    nomenclature_file: 'str'
      path to the nomenclature file (one available in besca/datasets)
    cnames: panda.DataFrame or list
        a dataframe with cluster to cell type attribution per distinct levels or a list of names
    reference_label: 'str'
        column in nomenclature file that corresponds to the names in cnames (short_dblabel if from make_annot function)
    new_label: 'str'
        column in nomenclature file that is used for the translation (e.g. 'dblabel' or 'EFO')
    new_level: int
        Keep the level if None or translate to a lower level e.g. 0, 1, 2, 3, 4

    returns
    -------
    panda.DataFrame containing the different levels of the nomenclature  indexed by cluster
    (based on cnames index)


    Example
    -------
    >>> import besca as bc
    >>> import pkg_resources
    >>> adata = bc.datasets.pbmc3k_processed()
    >>> new_cnames = bc.tl.sig.obtain_new_label(
    ...     nomenclature_file=pkg_resources.resource_filename('besca', 'datasets/nomenclature/CellTypes_v1.tsv'),
    ...     cnames=list(adata.obs['dblabel'].cat.categories),
    ...     reference_label='dblabel',
    ...     new_label = 'dblabel',
    ...     new_level = 2)

    """
    if type(cnames) == list:
        cnames = pd.DataFrame({"new_label": cnames})
        cnames = cnames.set_index(cnames["new_label"])
        cnames.index.name = None

    nomenclature = read_nomenclature(nomenclature_file)

    dblabel_dict = {}
    cat_dict = {}

    for mycol in list(cnames.columns):
        currentTerm = cnames[mycol]

        for cat in list(currentTerm):
            try:
                cat_nomenclature = nomenclature.loc[
                    nomenclature[reference_label] == cat, :
                ]
                # Initialize with dblabel, which is kept if is_level==new_level and if no better term can be found
                dblabel_dict[cat] = cat_nomenclature["dblabel"].iloc[0]
            except:
                raise Exception(
                    f"Error trying to reach {cat} in nomenclature. Please check if term exist"
                )

            if new_level is not None:
                is_level = int(cat_nomenclature["is_level"])
                # If is_level is greater than the new_level, try to find a better one between 0 and new_level, otherwise keep initial dblabel
                if is_level > int(new_level):
                    for l in range(0, int(new_level) + 1):
                        if not cat_nomenclature["l" + str(l)].isnull().values.any():
                            dblabel_dict[cat] = cat_nomenclature["l" + str(l)].iloc[0]
                # If is_level is smaller than the new_level, try to find a better one between is_level and new_level, otherwise keep initial dblabel
                elif is_level < int(new_level):
                    for l in range(is_level, int(new_level) + 1):
                        if not cat_nomenclature["l" + str(l)].isnull().values.any():
                            dblabel_dict[cat] = cat_nomenclature["l" + str(l)].iloc[0]

            # We retrieved the dblabel for each category. Let's translate to the new_label if needed:
            if new_label != "dblabel":
                cat_dict[cat] = nomenclature.loc[
                    nomenclature["dblabel"] == dblabel_dict[cat], new_label
                ].iloc[0]
            else:
                cat_dict = dblabel_dict

    new_cnames = []

    for mycol in list(cnames.columns):
        currentTerm = cnames[mycol]
        try:
            new_cnames.append(cat_dict[x] for x in list(currentTerm))
        except:
            raise Exception(
                f"Error trying to reach {currentTerm} in nomenclature. Please check if term exist"
            )

    new_cnames = pd.DataFrame(new_cnames).transpose()
    new_cnames.columns = cnames.columns
    new_cnames.index = cnames.index

    return new_cnames


def obtain_dblabel(nomenclature_file: str, cnames, reference_label="short_dblabel"):
    """Matches the cnames obtained by the make_annot function
    to the db label (standardized label).

    parameters
    ---------
    nomenclature_file: 'str'
      path to the nomenclature file (one available in besca/datasets)
    cnames: panda.DataFrame
        a dataframe with cluster to cell type attribution per distinct levels
    reference_label: 'str'
        column in nomenclature file that corresponds to the names in cnames (dblabel_short if from make_annot function)

    returns
    -------
    panda.DataFrame containing the different levels of the nomenclature  indexed by cluster
    (based on cnames index)
    """
    return obtain_new_label(
        nomenclature_file,
        cnames,
        reference_label=reference_label,
        new_label="dblabel",
        new_level=None,
    )


def read_nomenclature(nomenclature_file: str) -> pd.DataFrame:
    nomenclature = pd.read_csv(
        nomenclature_file, sep="\t", header=0, skiprows=range(1, 2)
    )
    return nomenclature


def match_label(
    vector_label: list,
    nomenclature_file: str,
    start_column: str = "dblabel",
    match_column: str = "short_dblabel",
) -> pd.DataFrame:
    """Return a table matching values in vector label.
    It is mean to translate labels , for example for plotting

      Parameters
      ----------
      vector_label: `list`
          initial list of values to match
      nomenclature_file: `str`
          location of the nomenclature_file
      start_column: `str` | default = dblabel
          column to start within the nomenclature
      match_column: `str` | default = short_dblabel
          column to match within the nomenclature

      Returns
      -------
      pd.DataFrame
          Pandas dataframe containing the matching labels. If matchings labels are not found, the value is kept as is.

      Example
      -------
      >>> import besca as bc
      >>> import scanpy as sc
      >>> import pkg_resources
      >>> adata = bc.datasets.pbmc3k_processed()
      >>> nomenclature_file=pkg_resources.resource_filename('besca', 'datasets/nomenclature/CellTypes_v1.tsv'),
      >>> nomenclature_file=''.join(nomenclature_file)
      >>> matching_v = bc.tl.sig.match_label(adata.obs.get( ["celltype3"]),  nomenclature_file)
      >>> adata.obs['shortlabel'] = adata.obs.get( "celltype3").map( dict(matching_v.values))
      >>> sc.pl.umap(adata, color=['shortlabel'])

    """
    nomenclature = read_nomenclature(nomenclature_file)
    if start_column not in nomenclature.columns:
        raise KeyError(
            f"{start_column} not found as one of the nomenclature file columns"
        )
    if match_column not in nomenclature.columns:
        raise KeyError(
            f"{match_column} not found as one of the nomenclature file columns"
        )

    matching_df = nomenclature[nomenclature[start_column].isin(vector_label)][
        [start_column, match_column]
    ].reset_index(drop=True)
    # check if values not in
    missed_values = np.isin(vector_label, nomenclature[start_column])
    if any(~missed_values):
        missed = np.array(vector_label)[~missed_values].tolist()
        print("Non found values in nomenclature files", set(missed))
        matching_df = pd.concat(
            [matching_df, pd.DataFrame({start_column: missed, match_column: missed})]
        )
    return matching_df
