import logging
import sys
import anndata
from numpy import ndarray, mean

import pytest

def frac_reads(adata):
    """Cacluate the fraction of reads being attributed to a specific gene.

    calculates the fraction of reads that are attributed to a specific gene over the entire dataset.
    This easily lets you identify genes which dominate all the reads in a dataset.

    parameters
    ----------
    adata: :class:`~anndata.AnnData`
        The annotated data matrix.
    threshold: `int` | default = 0
        Integer defining the value above which a gene is deemed as being expressed.

    returns
    -------
    None
        updates the:class:`~anndata.AnnData` object adding `frac_reads` to  `adata.var`.

    Example
    -------
    >>> import besca as bc
    >>> adata = bc.datasets.simulated_pbmc3k_raw()
    >>> bc.pp.frac_reads(adata)
    """

    if isinstance(adata, anndata.AnnData):

        if type(adata.X) == ndarray:
            adata_X = adata.X
        else:
            adata_X = adata.X.toarray()

        # calculate the fraction of reads that are attributed to each gene
        cum_reads_gene = sum(adata_X)
        frac_reads = cum_reads_gene / sum(cum_reads_gene)

        adata.var["frac_reads"] = frac_reads.tolist()

        return None

    else:
        sys.exit("please pass an AnnData object as data")


def frac_pos(adata, threshold=0):
    """Calculate the fraction of cells positive for expression of a gene.

    Calculates for each gene the fraction of cells that show a gene expression
    above the indicated threshold. The fraction of positive cells is added to
    adata.var

    parameters
    ----------
    adata: :class:`~anndata.AnnData`
        The annotated data matrix.
    threshold: `int` | default = 0
        Integer defining the value above which a gene is deemed as being expressed.

    returns
    -------
    None
        updates the:class:`~anndata.AnnData` object adding `frac_pos` to  `adata.var`.

    Example
    -------
    >>> import besca as bc
    >>> adata = bc.datasets.simulated_pbmc3k_raw()
    >>> bc.pp.frac_pos(adata)

    """
    if isinstance(adata, anndata.AnnData):

        if type(adata.X) == ndarray:
            adata_X = adata.X
        else:
            adata_X = adata.X.toarray()

        # calculate the percentage of cells that show an expression of the gene above the threshold
        fraction_pos = sum(adata_X > threshold) / adata.shape[0]

        # add calculated fraction positive to adata.var
        adata.var["frac_pos"] = fraction_pos.tolist()

        return None

    else:
        sys.exit("please pass an AnnData object as data")


def mean_expr(adata):
    """Calculate the mean expression of a gene.

    Calculates for each gene its mean expression in the adata object.

    parameters
    ----------
    adata: :class:`~anndata.AnnData`
        The annotated data matrix.

    returns
    -------
    None
        updates the:class:`~anndata.AnnData` object adding `mean` to  `adata.var`.

    Example
    -------
    >>> import besca as bc
    >>> adata = bc.datasets.simulated_pbmc3k_raw()
    >>> bc.pp.mean_expr(adata)

    """
    if isinstance(adata, anndata.AnnData):

        if type(adata.X) == ndarray:
            adata_X = adata.X
        else:
            adata_X = adata.X.toarray()

        # calculate the percentage of cells that show an expression of the gene above the threshold
        means = mean(adata_X, axis=0)

        # add calculated fraction positive to adata.var
        adata.var["mean"] = means.tolist()

        return None

    else:
        sys.exit("please pass an AnnData object as data")


def top_counts_genes(adata, top_n=10):
    """Give out the genes that contribute the largest fraction to the total UMI counts.

    Calculates the genes which contribute the largest fraction of UMI counts to the total
     UMI counts (i.e. the most abundant transcripts).

    parameters
    ----------
    adata: :class:`~anndata.AnnData`
        The annotated data matrix.
    top_n: `int` | default = 10
        Integer defining the number of entries to return. If None then the entire table is returned.

    returns
    -------
    DataFrame
        pandas DataFrame containing the top_n genes

    Example
    -------
    >>> import besca as bc
    >>> adata = bc.datasets.simulated_pbmc3k_raw()
    >>> genes = bc.pp.top_counts_genes(adata)

    """
    if isinstance(adata, anndata.AnnData):
        if "frac_reads" in adata.var.columns:

            # generate a table sorted by frac_pos (descending order)
            sorted = adata.var.sort_values(by="frac_reads", ascending=False)

            if top_n == None:
                return sorted
            else:
                return sorted.head(top_n)
        else:
            logging.info("calculating frac_reads")

            # calculate fraciton positive for a gene
            bdata = adata.copy()
            frac_reads(bdata)

            # generate a table sorted by frac_pos (descending order)
            sorted = bdata.var.sort_values(by="frac_reads", ascending=False)

            if top_n == None:
                return sorted
            else:
                return sorted.head(top_n)
    else:
        sys.exit("please pass an AnnData object as data")


def top_expressed_genes(adata, top_n=10):
    """Give out the genes most frequently expressed in cells.

    Calculates the genes most frequently expressed in cells and returns a table
    containing the top_n entries.

    parameters
    ----------
    adata: :class:`~anndata.AnnData`
        The annotated data matrix.
    top_n: `int` | default = 10
        Integer defining the number of entries to return. If None then the entire table is returned.

    returns
    -------
    DataFrame
        pandas DataFrame containing the top_n genes

    Example
    -------
    >>> import besca as bc
    >>> adata = bc.datasets.simulated_pbmc3k_raw()
    >>> genes = bc.pp.top_expressed_genes(adata)

    """
    if isinstance(adata, anndata.AnnData):
        if "frac_pos" in adata.var.columns:

            # generate a table sorted by frac_pos (descending order)
            sorted = adata.var.sort_values(by="frac_pos", ascending=False)

            if top_n == None:
                return sorted
            else:
                return sorted.head(top_n)
        else:
            logging.info("calculating frac_pos")

            # calculate fraciton positive for a gene
            bdata = adata.copy()
            frac_pos(bdata)

            # generate a table sorted by frac_pos (descending order)
            sorted = bdata.var.sort_values(by="frac_pos", ascending=False)

            if top_n == None:
                return sorted
            else:
                return sorted.head(top_n)
    else:
        sys.exit("please pass an AnnData object as data")
