# this module contains functions for cell type annotation based on signatures in python using scanpy
# import using the python version 1.3.2 at least is a prerequisit! needs to be checked in all functions

from numpy import array, log, where
from pandas import DataFrame, read_csv
from scanpy import AnnData
from scipy.stats import mannwhitneyu


def getset(df, signame_complete, threshold):
    """ Handles missing signatures aux function for make_anno
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

    if(bool(df.index.isin([signame_complete]).any())):
        return(set(df.columns[df.loc[signame_complete, :] > threshold]))
    else:
        return(set())


def gene_fr(f, Tcells):
    """ Returns the fraction of cells expressing a gene in distinct clusters
    Takes as an input a dataframe with fractions per clusters and a gene name
    Returns the fraction genes per cluster
    parameters
    ----------
    f: panda.DataFrame
      a dataframe of fraction genes per cell as outputted by besca
    Tcells: 'str'
      gene name

    returns
    -------
    panda.Series
        a Series of fractions per cluster and gene
    """
    k = f.loc[f['Description'].isin(Tcells), ].copy()
    k = k.iloc[:, 0:len(k.columns)]
    return(k.mean(axis=0, skipna=True))


def score_mw(f, mymarkers):
    """ Score Clusters based on a set of immune signatures to generate a df of pvals
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
    ids = set(f.iloc[:, 2:len(f.columns)])
    mypFrame = DataFrame(index=mymarkers.keys(), columns=list(ids))
    for key, value in mymarkers.items():
        for i in list(ids):
            mypFrame.loc[key, i] = -10*log(
                mannwhitneyu(x=f.loc[f['Description'].isin(value), :][i],
                             y=f[i], alternative='greater').pvalue
            )
    return(mypFrame)


def add_anno(adata, cnames, mycol, clusters='leiden'):
    """ Adds annotation generated with make_anno to a AnnData object
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

    notmatched = list(set(list(cnames.index))-set(newlab))
    if len(notmatched) > 0:
        errm = 'Please check that you have indexed the correct clustering. ' + \
            len(notmatched) + 'clusterIDs not found in ' + clusters
        return(errm)
    else:
        for i in list(set(newlab)):
            newlab = newlab.replace(i, cnames.loc[i, mycol]).copy()
        return(newlab)


def read_annotconfig(configfile):
    """ Reads the configuration file for signature-based hierarhical annotation.

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
    sigconfig = read_csv(configfile, sep='\t', index_col=0)
    # Reorder with the specified order. Place better signatures first, only first match will be kept
    sigconfig = sigconfig.sort_values('Order')
    # Consider up to 7 levels
    nochild = list(set(sigconfig.index)-set(sigconfig['Parent']))
    levs = []
    levs.append(list(sigconfig.loc[sigconfig['Parent'] == 'None', ].index))
    for i in range(0, 7):
        levs.append(
            list(sigconfig.loc[sigconfig['Parent'].isin(levs[i]), ].index))

    levsk = []
    for lev in levs:
        if len(lev) > 0:
            levsk.append(lev)

    return(sigconfig, levsk)



def export_annotconfig(sigconfig, levsk, resultsFolder, analysisName = ''):
    """ Export the configuration defined in sigconfig and levsk
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
    it_order = 0 # Iterative order.
    config_t = sigconfig.copy()
    get_order = {}
    for child in levsk :
      for element in child :
              get_order[element] = it_order
              it_order +=1
    config_t['Order'] = config_t.index.map(get_order)
    config_t.to_csv( resultsFolder + analysisName + '_config.tsv' ,  sep = '\t' )






def make_anno(df, sigscores, sigconfig, levsk, lab='celltype', toexclude=[]):
    """ Annotate cell types
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
    for x in levsk:
        tmp = []
        for y in x:
            toinc = set(sigscores.keys())-set(toexclude)
            if y in list(toinc):
                tmp.append(y)
        levskk.append(tmp)
    levsk = levskk.copy()

    sigscoresk = {}
    for x in sigscores.keys():
        if not x in toexclude:
            sigscoresk[x] = sigscores[x]
    sigscores = sigscoresk.copy()

    # First part, get cluster identities
    myclust = list(df.columns)
    annol = {}
    for ii in myclust:
        annol[ii] = []
        for lev in levsk[0]:
            if ii in sigscores[lev]:
                annol[ii] = lev
                break
        if (len(levsk) > 2):
            for j in range(len(levsk)-1):
                # jj=j+1
                if len(annol[ii]) == 0:
                    break
                if len(str.split(annol[ii], ".")) < (j+1):
                    break
                sublev = list(
                    sigconfig.loc[sigconfig["Parent"] == str.split(annol[ii], ".")[j], :].index)
                sublevk = []
                for x in sublev:
                    if x in set(list(sigscores.keys())):
                        sublevk.append(x)
                sublev = sublevk.copy()
                for lev in sublev:
                    if ii in sigscores[lev]:
                        annol[ii] = annol[ii]+'.'+lev
                        break

    # Second part, transform to pandas
    cnames = []
    for key, val in annol.items():
        if len(val) == 0:
            tmp = ['Cell']
        else:
            tmp = str.split(val, ".")
        for i in range(len(levsk)-len(tmp)):
            tmp.append(tmp[len(tmp)-1])
        cnames.append(tmp)
    cnames = DataFrame(cnames)
    cnames.index = annol.keys()
    cnames.columns = [lab+str(x) for x in list(cnames.columns)]

    return(cnames)


def match_cluster(adata, obsquery, obsqueryval, obsref='leiden', cutoff=0.5):
    """ Matches categories from adata.obs to each other. For a query category specified in obsquery
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
    myleiden = list(
        set(adata[adata.obs[obsquery].isin([obsqueryval])].obs[obsref]))
    for i in myleiden:
        mysub = adata[adata.obs[obsref].isin([i])].copy()
        myc.append(len(mysub[mysub.obs[obsquery] == obsqueryval])/len(mysub))
    return(list(array(myleiden)[array(myc) > cutoff]))


def obtain_dblabel(nomenclature_file, cnames):
    """ Matches the cnames obtained by the make_annot function
    to the db label (standardized label).

    parameters
    ---------
    nomenclature_file: 'str'
      path to the nomenclature file (one available in besca/datasets)
    cnames: panda.DataFrame
        a dataframe with cluster to cell type attribution per distinct levels

    returns
    -------
    panda.DataFrame containing the different levels of the nomenclature  indexed by cluster
    (based on cnames index)
    """
    nomenclature = read_csv(nomenclature_file, sep='\t',
                            header=0, skiprows=range(1, 2))
    cnamesDBlabel = []
    for mycol in list(cnames.columns):
        currentTerm = cnames[mycol]
        try:
            cnamesDBlabel.append([list(nomenclature.loc[nomenclature['short_dblabel'] == x, 'dblabel'])[
                                 0] for x in list(currentTerm)])
        except:
            raise Exception(
                f'Error trying to reach {currentTerm} in nomenclature. Please check if term exist')
    cnamesDBlabel = DataFrame(cnamesDBlabel).transpose()
    cnamesDBlabel.columns = cnames.columns
    cnamesDBlabel.index = cnames.index
    return(cnamesDBlabel)
