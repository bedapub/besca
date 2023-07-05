import besca as bc
import re
import pandas as pd
import pytest

def filter_perturb(adata, col = "n_sgRNAs"):
    """
    Filters cells based on number of perturbations.

    parameters

    ----------
    adata: `AnnData`
        AnnData object containing data that is to be filtered.
    col: `str`
        The adata.obs collumn that contains the number of sgRNAs that have perturbed the cell

    returns
    -------

    adata: Filtered adata, with cells having only 1 perturbation

    Examples
    --------

    >>> pytest.skip('Test will be skipped, because the dataset is not available on zenodo anymore.')
    >>> import besca as bc
    >>> adata = bc.datasets.crispr_10x_unfiltered()
    >>> adata = bc.pp.filter_perturb(adata, col = "n_sgRNAs")
    """


    filters = adata.obs[col] > 1

    #Make sure cells have perturbations of more than 1
    if True in filters.value_counts():
        adata = bc.subset_adata(adata, filter_criteria= ~filters, axis = 0, raw=False)
        print("removed {} with more than one perturbation".format(filters.value_counts()[True]))

    #Get the cells with different than 0 perturbations
    filters = adata.obs[col] != 0
    adata = bc.subset_adata(adata, filter_criteria= filters, axis = 0, raw=False)

    print("removed {} with no perturbation".format(filters.value_counts()[False]))

    return adata

def extract_target(adata, col = "sgRNAs", col_target = "Target", col_id = "samples__sample_id"):
    """
    Create a new variable in obs that contains the gene-target per perturbation

    parameters

    ----------
    adata: `AnnData`
        AnnData object containing the sgRNAs per cell.
    col: `str`
        The adata.obs collumn that contains the sgRNAs that have perturbed the cell
    col_target: `str`
        Prefered name of gene-target column
    col_id: `str`
        Name of the column that contains the samples id

    returns
    -------

    adata: Annotated data with an extra variable in obs

    Examples
    --------

    >>> pytest.skip('Test will be skipped, because the dataset is not available on zenodo anymore.')
    >>> import besca as bc
    >>> import pandas as pd
    >>> import re
    >>> adata = bc.datasets.crispr_10x_filtered()
    >>> adata = bc.pp.extract_target(adata)
    """
    adata.obs[col_target] = "Target"
    def set_target(row):
        gRNA = row[col]
        if "No-Guide" == gRNA:
            return gRNA
        elif "Target" in gRNA:
            return "Control"
        if "," in gRNA:
            gRNA = gRNA.split(",")[0]
        gRNA = re.split(r'[_|-]',gRNA)[0]
        return gRNA

    #For each row based on the sgRNA, set the name of the targeted gene
    adata.obs[col_target] = adata.obs.apply(set_target, axis = 1)
    perturbations_per_sample = pd.crosstab(index=adata.obs[col_id], columns=adata.obs[col_target])

    #Remove perturbations with low level of infection
    cells_to_kick = {}
    for samples,rows in perturbations_per_sample.iterrows():
        cells_to_kick[samples] = []
        for name, value in rows.iteritems():
            if value <=5 and value >0:
                cells_to_kick[samples].append(name)

    def remove_cells_per_sample(row):
        if row[col_target] in cells_to_kick[row[col_id]]:
            return False
        else:
            return True
    my_filter = adata.obs.apply(remove_cells_per_sample, axis = 1)
    adata = bc.subset_adata(adata, filter_criteria= my_filter, axis = 0, raw=False)

    return adata
