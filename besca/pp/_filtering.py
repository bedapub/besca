import sys
import anndata
from pandas import read_csv
from numpy import sum
from besca.datasets._mito import get_mito_genes

# import helper functions from besca
from besca._helper import get_raw


def filter(
    adata,
    max_genes=None,
    min_genes=None,
    max_counts=None,
    min_counts=None,
    min_cells=None,
    max_mito=None,
    annotation_type=None,
    species="human",
):
    """Filter cell outliers based on counts, numbers of genes expressed, number of cells expressing a gene and mitochondrial gene content.

    Filtering is performed iteratively in the order: max_genes, min_genes, max_counts, min_counts,
    min_cells, max_mito.

    The Thresholds are defined as follows:
    max_genes >= n_genes
    min_genes <= n_genes
    max_counts >= n_counts
    min_counts <= n_counts
    min_cells <= n_cells
    max_mito > percent mito

    parameters
    ----------
    adata: :class:`~anndata.AnnData`
        The annotated data matrix.
    max_genes: `int` | default = None
        integer value specifying the threshold for the maximum number of genes that a cell needs to express
    min_genes: `int` | default = None
        integer value specifying the threshold for the minimum number of genes that a cell needs to express
    max_counts: `int` | default = None
        integer value specifying the maximum number of counts that a cell needs to contain
    min_counts: `int` | default = None
        integer value specifying the minimum number of counts that a cell needs to contain
    min_cells: `int` | default = None
        integer value specifying the minimum number of cells that need to express a gene for it to be included
    max_mito: `float` | default = None
        decimalvalue specifying the threshold for the maximum percentage of mitochondrial genes in a cell
    annotation_type: `SYMBOL` or `ENSEMBLE` or `None` | default = None
        string identifying the type of gene ids contained in adata.var_names, necessary for identifying
        mitochondrial genes in case percent_mito is not already included in adata.obs

    returns
    -------
    AnnData
        filters the :class:`~anndata.AnnData` object and potentially adds either `n_genes` or `n_counts` to `adata.obs`.

    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.pbmc3k_raw()
    >>> adata.size
    >>> adata = bc.pp.filter(adata, max_counts=6500, max_genes=1900, max_mito=0.05, min_genes=600, min_counts=600, min_cells=2)
    >>> adata.size


    """
    if isinstance(adata, anndata.AnnData):

        # get original number of cells
        ncells = adata.shape[0]
        ngenes = adata.shape[1]
        print(
            "started with ",
            str(ncells),
            " total cells and ",
            str(ngenes),
            " total genes",
        )

        # calculate values if necessary
        if max_counts is not None or min_counts is not None:
            if adata.obs.get("n_counts") is None:
                adata.obs["n_counts"] = adata.X.sum(axis=1)
        if min_genes is not None or max_genes is not None:
            if adata.obs.get("n_genes") is None:
                adata.obs["n_genes"] = sum(adata.X > 0, axis=1)
        if min_cells is not None:
            if adata.var.get("n_cells") is None:
                adata.var["n_cells"] = sum(adata.X > 0, axis=0).T

        if max_mito is not None:
            if adata.obs.get("percent_mito") is None:
                mito_list = get_mito_genes(species, annotation_type)
                mito_genes = [
                    name for name in adata.var_names if name in mito_list
                ]  # ensembl
                # for each cell compute fraction of counts in mito genes vs. all genes
                n_counts = sum(adata.X, axis=1).A1
                n_counts[n_counts == 0] = float("inf")
                adata.obs["percent_mito"] = (
                    sum(adata[:, mito_genes].X, axis=1).A1 / n_counts
                )

        # perform actual filtering for all given parameters
        if max_genes is not None:
            curr_cells = adata.shape[0]
            adata = adata[adata.obs.get("n_genes") <= max_genes, :].copy()
            new_cells = adata.shape[0]
            print(
                "removed",
                str(curr_cells - new_cells),
                "cells that expressed more than",
                str(max_genes),
                "genes",
            )
        if min_genes is not None:
            curr_cells = adata.shape[0]
            adata = adata[adata.obs.get("n_genes") >= min_genes, :].copy()
            new_cells = adata.shape[0]
            print(
                "removed",
                str(curr_cells - new_cells),
                "cells that did not express at least",
                str(min_genes),
                " genes",
            )
        if max_counts is not None:
            curr_cells = adata.shape[0]
            adata = adata[adata.obs.get("n_counts") <= max_counts, :].copy()
            new_cells = adata.shape[0]
            print(
                "removed",
                str(curr_cells - new_cells),
                "cells that had more than",
                str(max_counts),
                " counts",
            )

        if min_counts is not None:
            curr_cells = adata.shape[0]
            adata = adata[adata.obs.get("n_counts") >= min_counts, :].copy()
            new_cells = adata.shape[0]
            print(
                "removed",
                str(curr_cells - new_cells),
                "cells that did not have at least",
                str(min_counts),
                "counts",
            )

        if min_cells is not None:
            curr_genes = adata.shape[1]
            adata = adata[:, adata.var.get("n_cells") >= min_cells].copy()
            new_genes = adata.shape[1]
            print(
                "removed",
                str(curr_genes - new_genes),
                "genes that were not expressed in at least",
                str(min_cells),
                "cells",
            )

        if max_mito is not None:
            curr_cells = adata.shape[0]
            adata = adata[adata.obs.get("percent_mito") < max_mito, :].copy()
            new_cells = adata.shape[0]
            print(
                "removed ",
                str(curr_cells - new_cells),
                " cells that expressed ",
                str(max_mito * 100),
                "percent mitochondrial genes or more",
            )

        ncells_final = adata.shape[0]
        ngenes_final = adata.shape[1]
        print(
            "finished with",
            str(ncells_final),
            " total cells and",
            str(ngenes_final),
            "total genes",
        )

        return adata

    else:
        print("please pass an AnnData object as data")
        sys.exit(1)


def filter_gene_list(adata, filepath, use_raw=True, use_genes="SYMBOL"):
    """Function to remove all genes specified in a gene list read from file.

    This function removes all genes in the dataset from a given list of genes (for example mito genes or ribosomal genes).
    Note that the input file consists of two columns (ENSEMBL gene id and gene symbol).

    Parameters
    ----------
    adata: `AnnData`
        AnnData object
    filepath: `str`
        File path as string point to gene list
    use_raw: `bool` | default = True
        Boolian indicator if the function should return the filtered raw object or not.
    use_genes: `SYMBOL` or `ENSEMBL` | default = SYMBOL
        String defining whether ENSEMBL id's or gene symbols are used in the adata.var_names (defines which column of input gene list is read)

    Returns
    -------
    AnnData
        Returns the filtered AnnData object.

    Example
    -------
    >>> import besca as bc
    >>> import os
    >>> adata = bc.datasets.pbmc_raw()
    >>> adata.n_vars
    >>> # TO COMPLETE reference_mito_file = 'path2file.tsv'
    >>> # adata = bc.pp.filter_gene_list(adata, file_path=reference_mito_file, use_genes='SYMBOL')
    >>> # adata.n_vars

    """
    # generate copy of adata object
    if use_raw:
        if adata.raw is None:
            print(
                "WARNING: adata does not contain .raw filtering on regular adata object"
            )
            adata = adata.copy()
        else:
            adata = get_raw(adata)
    else:
        adata = adata.copy()

    if use_genes == "SYMBOL":
        gene_list = list(
            read_csv(filepath, header=None, sep="\t")[1]
        )  # ENS_GENE_ID  GENE_SYMBOL (2 cols)
    elif use_genes == "ENSEMBL":
        gene_list = list(
            read_csv(filepath, header=None, sep="\t")[0]
        )  # ENS_GENE_ID  GENE_SYMBOL (2 cols)
    else:
        sys.exit("Please supply either SYMBOL or ENSEMBL ids")
    genes = [i for i in adata.var_names if i not in gene_list]
    # remove all the genes from the list
    if len(genes) > 0:
        gene_indexes = [adata.var_names.tolist().index(x) for x in genes]
        adata = adata[:, gene_indexes]
    else:
        print(
            "None of the genes from input list found in data set. Please ensure you have correctly specified use_genes to match the type of genes saved in adata.var_names."
        )
    return adata
