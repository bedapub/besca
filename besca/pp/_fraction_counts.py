import sys
from pandas import read_csv
from numpy import sum, any
import warnings
from besca.datasets._mito import get_mito_genes


def fraction_counts(
    adata, species="human", name="percent_mito", use_genes="SYMBOL", specific_file=None
):
    """Function to calculate fraction of counts per cell from a gene list.
    This function calculates the fraction of counts per cell for
    a list of genes (for example mito genes) if no specific file is given.
    Note that the input file consists of two columns
    (ENSEMBL gene id and gene symbol) tab
    separated

    Parameters
    ----------
    adata: `AnnData`
        AnnData object
    species: `str` | default = human
        species for mitochondrial content evaluation
    name: `str` | default = percent_mito
        String identifying the column name to which the results
        should be written to in the AnnData.obs object
    use_genes: `SYMBOL` or `ENSEMBL` | default = SYMBOL
        String defining whether ENSEMBL id's or gene symbols are used in the
        adata.var_names (defines which column of input gene list is read)
    specific_file: `str` | default None.
        if indicated, the file will be used to extract the gene list
    Returns
    -------
    None
        Returns None but updates adata with new column named 'name'
        containing calculated fraction of counts.

    Example
    -------

    >>> import besca as bc
    >>> import os
    >>> adata = bc.datasets.pbmc_raw()
    >>> adata.obs.head(5)
    >>> bc.pp.fraction_counts(adata,  'human', use_genes='SYMBOL')
    >>> adata.obs.head(5)
    """
    if specific_file is None:
        gene_list = get_mito_genes(species, use_genes)
    else:
        # ENS_GENE_ID  GENE_SYMBOL (2 cols)
        if use_genes == "SYMBOL":
            gene_list = list(read_csv(specific_file, header=None, sep="\t")[1])
        elif use_genes == "ENSEMBL":
            gene_list = list(read_csv(specific_file, header=None, sep="\t")[0])
    genes = [i for i in adata.var_names if i in gene_list]
    # for each cell compute fraction of counts in gene_list vs. all genes
    # axis=1 --> sum over rows
    if len(genes) > 0:
        n_counts = sum(adata.X, axis=1).A1
        if any(n_counts == 0):
            warnings.warn(
                "Some of the cells contain no counts. \
                           Do not forget to remove 'empty' cells from data."
            )
            n_counts[n_counts == 0] = float("inf")
        adata.obs[name] = sum(adata[:, genes].X, axis=1).A1 / n_counts
    else:
        adata.obs[name] = 0.0
        print(
            "None of the genes from input list found in data set. \
               Please ensure you have correctly specified use_genes to match \
               the type of genes saved in adata.var_names."
        )
    return None
