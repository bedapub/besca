import math
import besca as bc
import pandas as pd
import pytest
from scipy.spatial import distance

def execute_de_sgRNA(adata, mypairs, column_focus, myfc, mypval):
    """
    Performs DGEX, based on pairs of values that exist in the "column_focus"

    parameters

    ----------
    adata: `AnnData`
        AnnData object containing the sgRNAs per cell.
    mypairs: list`
        pairs of labels to do the DGEX
    column_focus: `str`
       The column that the pair values exist in
    myfc: `float`
        Fold change threshold
    mypval: `float`
        P-value threshold

    returns
    -------

    DGEX as many times as the pairs

    Examples
    --------

    >>> pytest.skip('Test will be skipped, because the dataset is not available on zenodo anymore.')
    >>> import besca as bc
    >>> import math
    >>> import pandas as pd
    >>> adata = bc.datasets.crispr_10x_filtered()
    >>> dgex = bc.tl.crispr.execute_de_sgRNA(adata, [('RAB1A','Control')], 'Target', 0, 1)
    """
    #DGEX for each of the pairs of values
    DElist={}
    for i in mypairs:
        DElist[i[0]+'-'+i[1]]=bc.tl.dge.get_de(adata[adata.obs[column_focus].isin(i)],column_focus,demethod='wilcoxon',topnr=20000, logfc=myfc,padj=mypval)
    dge_data =[]

    #Prepare to store DGEX data in dataframe
    for i in mypairs:
        dge_list_c = DElist[i[0]+'-'+i[1]][i[1]]
        dge_list_t = DElist[i[0]+'-'+i[1]][i[0]]
        if not dge_list_c[dge_list_c.Name == i[0]].empty:
            dge_data.append(dge_list_c[dge_list_c.Name == i[0]])
        elif not dge_list_t[dge_list_t.Name == i[0]].empty:
            dge_list_t.loc[dge_list_t.Name == i[0],"Log2FC"] = -dge_list_t[dge_list_t.Name == i[0]].Log2FC
            dge_data.append(dge_list_t[dge_list_t.Name == i[0]])
        else:
            data = [[i[0],0,0,0]]
            dge_data.append(pd.DataFrame(data, columns=['Name','Score','Log2FC','P.adj']))
    dge_data = pd.concat(dge_data).reset_index(drop=True).sort_values(by=['Log2FC'], ascending=False)

    # Turn padj to log padj and log fold change in respect to Control
    dge_data = dge_data.apply(_prepare,axis=1)
    return dge_data

def _prepare(row):
    row['Log2FC'] = -row['Log2FC']
    if row['P.adj'] != 0:
        row['P.adj'] = float(math.log10(row['P.adj']))
        row['P.adj'] = -row['P.adj']
    return row

def find_distances(adata, experiment = "CROPseq"):
    """
    Distance matrix between cells using the PCA results

    parameters

    ----------
    adata: `AnnData`
        AnnData object containing the sgRNAs per cell after conducting PCA.
    experiment: `str`
        Type of Crispr Experiment
    returns
    -------
    The Distance Matrix using Euclidean  metric

    Examples
    --------

    >>> pytest.skip('Test will be skipped, because the dataset is not available on zenodo anymore.')
    >>> import besca as bc
    >>> import math
    >>> import pandas as pd
    >>> from scipy.spatial import distance
    >>> adata = bc.datasets.crispr_10x_filtered()
    >>> adata_subset = bc.tl.crispr.find_distances(adata, experiment = '10xChromium')
    """

    # Data containing the identifiers
    if all(x in adata.obs.columns for x in ['sgRNAs','Target','leiden']):

        if experiment == "CROPseq":
            df=adata.obs[['group', 'sgRNAs','Target','leiden']].copy()
        else:
            df=adata.obs[['sgRNAs','Target','leiden']].copy()

        #Compute distance matrixes and store in a dataframe
        Y = distance.cdist(adata.obsm['X_pca'][:,:2], adata.obsm['X_pca'][:,:2], 'euclidean')
        Euc_dist = pd.DataFrame(Y, columns=df.index, index=df.index)

        return Euc_dist
    else:
        print("Necessary columns not present. Necessary columns are ['sgRNAs','Target','leiden']")
        return None
