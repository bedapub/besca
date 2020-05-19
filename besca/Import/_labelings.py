from pandas import read_csv
from os.path import join as pathjoin

def add_cell_labeling(adata,
                      filepath,
                      label = 'celltype'):
    """ add a labeling written out in the FAIR formating to adata.obs

    A laveling contained in the FAIR compliant cell2labels.tsv is added to adata.obs. The string supplied
    in the parameter label is used to name the column in adata.obs that contains the imported labeling.

    All cells that are not labeled in the cell2labels.tsv will be annotated with 'not labeled'.
    
    parameters
    ----------
    adata: `AnnData`
        the AnnData object whose obs should be updated
    filepath: `str` 
        filepath to the cell2labels.tsv that is to be appended to adata.obs
    label: `str` | default = 'celltype'
        string indicating the label that is to be added to the annotation that is being imported
    
    returns
    -------
    None
        updates the supplied AnnData object

    """

    labeling = read_csv(pathjoin(filepath, 'cell2labels.tsv'), sep='\t', index_col=0)
    labeling.rename(columns={'LABEL':label}, inplace=True)

    adata.obs[label] = 'not labeled'
    adata.obs.update(labeling)

    return(None)
