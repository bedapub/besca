from scanpy import AnnData
from scanpy.external.pp import mnn_correct
from pandas import DataFrame
import scipy

def batch_correct(adata, batch_to_correct):
    """function to perform batch correction

    mnnpy batch correction with the batch specified by batch_to_correct using all the
    genes contained in adata. Before running this function the highly variable genes should be selected.

    This function will return a corrected AnnData object without any .raw! please use the function postprocess_mnnpy to generate
    a complete AnnData object which contains the corrected values and the original raw. Because of this pleae ensure that you do
    not save the output of this function to the variable adata!

    parameters
    ----------
    adata: AnnData
        AnnData objec that is to be batch corrected
    batch_to_correct: 'str'
        string identifiyng the column which contains the batches that are to be corrected

    returns:
    AnnData
        the corrected AnnData object (raw values have been truncated to the highly variable genes)
    """
    batches = list(set(adata.obs[batch_to_correct]))
    batchNum = len(batches)
    batchData = []
    for i in range(batchNum):
        batchData.append(adata[adata.obs [batch_to_correct] == batches [i]].copy())
    corrected = mnn_correct(*batchData, var_index = None, batch_key = batch_to_correct, do_concatenate = True, svd_mode='svd', save_raw = False)
    (corrected [0]).var = adata.var.copy()
    # otherwise, the symbols and n_cells columns will be identical copies for each batch, with modified names
    return(corrected[0])


def postprocess_mnnpy(adata, bdata):
    """ postprocessing to generate a newly functional AnnData object

    After running mnnpy_mnncorrect we obtain ann AnnData object bdata. Since mnn_correct automatically
    truncates all the genes contained in .raw to contain only the highly variable genes this function
    creates a new AnnData object that contains .X from bdata but .raw from AnnData (which still contains all the
    genes, not only the highly variable ones).

    Before creation of the new AnnData object the matrices are sorted according to cellbarcode so
    that we ensure the labelings are correct.

    parameters
    ----------

    adata:
        the uncorrected AnnData object
    bdata:
        the batch correted AnnData object

    returns
    -------
    AnnData
        AnnData object with adata.X containing the corrected values and .raw all of the original values

    """
    corrected_matrix = DataFrame(data = bdata.X, index = bdata.obs_names.tolist(), columns = bdata.var_names.tolist())
    corrected_matrix.sort_index(inplace=True)

    new_adata = AnnData(corrected_matrix.values)
    new_adata.obs = bdata.obs.sort_index()
    new_adata.var_names = bdata.var_names
    new_adata.obs_names = bdata.obs_names.sort_values()
    new_adata.var = bdata.var

    #need to sort raw object to match the batch corrected order
    raw_matrix = DataFrame(data=(adata.raw.X.todense() if scipy.sparse.issparse(adata.raw.X) else adata.raw.X), index=adata.obs_names.tolist(), columns=adata.raw.var_names.tolist())
    raw_matrix.sort_index(inplace=True)

    #recreate raw
    raw = AnnData(raw_matrix.values)
    raw.var_names = adata.raw.var_names
    raw.obs_names = adata.obs_names.sort_values()
    raw.var = adata.raw.var

    #add raw back in
    new_adata.raw = raw

    #ensure that indices are preserved
    adata.obs_names = adata.obs.CELL
    adata.obs.index = adata.obs.CELL

    return(new_adata)
