import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import anndata
from besca.datasets._mito import get_mito_genes

def kp_genes(adata, threshold=0, min_genes=100, ax=None, figsize=None):
    """visualize the minimum gene per cell threshold.

    Generates a "knee-plot" to visualize the chosen threshold for the minimum number of genes that a cell expresses to validate if this is a good threshold.

    Parameters
    ----------
    adata: :class:`~anndata.AnnData`
        The annotated data matrix.
    threshold: `int` | default = 0
        integer value that defines the minimum expression threshold for a gene to be defined as expressed. Default value is 0.
    min_genes: `int` | default = 100
        visualize the chosen minimum gene threshold (default is set to 100)
    ax: `axes` | default = None
        pass the axes class to which your figure should be added
    figsize: (width, height) or None | default = None
        optional parameter to define the figure size of the plot that is to be generated

    Returns
    -------
    None
        Figure is displayed

    Example
    -------

    Generate a "knee-plot" for a minimum of 600 expressed genes per cell for an example dataset.
    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> adata = bc.datasets.simulated_pbmc3k_raw()
    >>> min_genes = 600
    >>> fig, ax1 = plt.subplots(1)
    >>> bc.pl.kp_genes(adata, min_genes = min_genes, ax = ax1)

    .. plot::
        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> adata = bc.datasets.simulated_pbmc3k_raw()
        >>> min_genes = 600
        >>> fig, ax1 = plt.subplots(1)
        >>> bc.pl.kp_genes(adata, min_genes = min_genes, ax = ax1)

    """
    nbr_cells = adata.shape[0]
    ax = ax or plt.gca()

    ax.plot(-np.sort(-np.array(np.sum(adata.X > threshold, axis=1).T)[0]))
    ax.set_yscale("log", base=10)
    ax.set_xscale("log", base=10)
    ax.set_ylabel("Number of expressed genes")
    ax.set_xlabel("Cells")
    ax.set_title("Expressed genes [count > " + str(threshold) + "]")
    if figsize is not None:
        plt.figure(figsize=figsize)
    ax.hlines(min_genes, xmin=0, xmax=nbr_cells * 1.3, color="red", linestyles="dotted")

    return None


def kp_counts(adata, min_counts=200, ax=None, figsize=None):
    """visualize the minimum UMI counts per cell threshold.

    this function generates a knee-plot visualizing a given min_counts cutoff when given an adata object

    Parameters
    ----------
    adata: :class:`~anndata.AnnData`
        The annotated data matrix.
    min_counts: `int` | default = 200
        visualize the chosen minimum UMI counts threshold
    ax: `axes` | default = None
        pass the axes class to which your figure should be added
    figsize: (width, height) or None | default = None
        optional parameter to define the figure size of the plot that is to be generated


    Returns
    -------
    None
        figure is displayed


    Example
    -------

    Generate a "knee-plot" for a minimum of 600 UMI counts per cell for an example dataset.
    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> adata = bc.datasets.simulated_pbmc3k_raw()
    >>> min_counts = 600
    >>> fig, ax1 = plt.subplots(1)
    >>> bc.pl.kp_counts(adata, min_counts = min_counts, ax = ax1)

    .. plot::
        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> adata = bc.datasets.simulated_pbmc3k_raw()
        >>> min_counts = 600
        >>> fig, ax1 = plt.subplots(1)
        >>> bc.pl.kp_counts(adata, min_counts = min_counts, ax = ax1)

    """
    nbr_cells = adata.shape[0]
    ax = ax or plt.gca()

    ax.plot(-np.sort(-np.array(adata.X.sum(axis=1).T)[0]))
    ax.set_yscale("log", base=10)
    ax.set_xscale("log", base=10)
    ax.set_ylabel("Number of UMI counts")
    ax.set_xlabel("Cells")
    ax.set_title("Counts per cell")
    if figsize is not None:
        plt.figure(figsize=figsize)
    ax.hlines(
        min_counts, xmin=0, xmax=nbr_cells * 1.3, color="red", linestyles="dotted"
    )

    return None


def kp_cells(adata, threshold=0, min_cells=2, ax=None, figsize=None):

    """visualize the minimum number of cells expressing a gene threshold.

    This function generates a "knee-plot" visualizing a given min_cells cutoff when given an adata object.
    Threshold sets the value above which a gene is defined as being expressed. All of the genes to the right
    of the red vertical line would no longer be included in the dataset if you choose to filter with the shown
    parameter.

    Parameters
    ----------
    adata: :class:`~anndata.AnnData`
        The annotated data matrix.
    threshold: `int` | default = 0
        integer value that defines the minimum expression threshold above which a gene to be defined as expressed. Default value is 0.
    min_cells: `int` | default = 2
        visualize the chosen minimum number of cells that need to express a gene threshold
    ax: `axes` | default = None
        pass the axes class to which your figure should be added
    figsize: (width, height) or None | default = None
        optional parameter to define the figure size of the plot that is to be generated

    Returns
    -------
    None
        Figure is displayed

    Example
    -------

    Generates a "knee-plot" for a the minimum number of cells expressing a gene in a given dataset.
    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> adata = bc.datasets.simulated_pbmc3k_raw()
    >>> min_cells = 2
    >>> fig, ax1 = plt.subplots(1)
    >>> bc.pl.kp_cells(adata, min_cells = min_cells, ax = ax1)

    .. plot::
        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> adata = bc.datasets.simulated_pbmc3k_raw()
        >>> min_cells= 2
        >>> fig, ax1 = plt.subplots(1)
        >>> bc.pl.kp_cells(adata, min_cells = min_cells, ax = ax1)


    """
    nbr_cells = adata.shape[0]

    ax = ax or plt.gca()

    ax.plot(-np.sort(-np.array(np.sum(adata.X > threshold, axis=0))[0]))
    ax.set_yscale("log", base=10)
    ax.set_xscale("log", base=10)
    ax.set_xlabel("genes")
    ax.set_ylabel("number of cells expressing a gene")
    ax.set_title("Cells that express a gene [count > " + str(threshold) + "]")

    cutoff = np.sum(
        -np.sort(-np.array(np.sum(adata.X > threshold, axis=0))) > min_cells
    )
    if figsize is not None:
        plt.figure(figsize=figsize)
    ax.vlines(cutoff, 0, nbr_cells * 1.05, color="red", linestyle="dotted")

    return None


def max_counts(adata, max_counts=10000, ax=None, figsize=None):
    """visualize maximum UMI counts per cell threshold.

    this function generates a knee-plot visualizing a given min_cells cutoff when given an adata object

    Parameters
    ----------
    adata: :class:`~anndata.AnnData`
        The annotated data matrix.
    max_counts: `int` | default = 10000
        visualize the chosen maximum UMI counts cutoff
    ax: `axes` | default = None
        pass the axes class to which your figure should be added
    figsize: (width, height) or None | default = None
        optional parameter to define the figure size of the plot that is to be generated

    Returns
    -------
    None
        Figure is displayed

    Example
    -------

    Generates a "knee-plot" for a maximum UMI count of 6500 for the example dataset.
    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> adata = bc.datasets.simulated_pbmc3k_raw()
    >>> max_counts = 6500
    >>> fig, ax1 = plt.subplots(1)
    >>> bc.pl.max_counts(adata, max_counts = max_counts, ax = ax1)

    .. plot::
        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> adata = bc.datasets.simulated_pbmc3k_raw()
        >>> max_counts = 6500
        >>> fig, ax1 = plt.subplots(1)
        >>> bc.pl.max_counts(adata, max_counts = max_counts, ax = ax1)

    """
    if adata.obs.get("n_counts") is None:
        adata.obs["n_counts"] = adata.X.sum(axis=1)
    if adata.obs.get("n_genes") is None:
        adata.obs["n_genes"] = np.sum(adata.X > 0, axis=1)

    x = adata.obs.get("n_counts").tolist()
    y = adata.obs.get("n_genes").tolist()

    data_plot = pd.DataFrame({"n_counts": x, "n_genes": y})
    min_genes = 100

    ax = ax or plt.gca()

    ax.scatter(
        data=data_plot, x="n_counts", y="n_genes", alpha=0.4, s=0.5, rasterized=True
    )
    if figsize is not None:
        plt.figure(figsize=figsize)
    ax.vlines(max_counts, min_genes, max(y), color="red", linestyle="dotted")
    ax.set_xlabel("n_counts")
    ax.set_ylabel("n_genes")
    ax.set_title("max counts cutoff")

    return None


def max_genes(adata, max_genes=5000, ax=None, figsize=None):
    """visualize maximum number of genes per cell threshold.

    this function generates a knee-plot visualizing a given min_cells cutoff when given an adata object

    Parameters
    ----------
    adata: :class:`~anndata.AnnData`
        The annotated data matrix.
    max_genes: `int` | default = 5000
        visualize the chosen maximum gene number cutoff
    ax: `axes` | default = None
        pass the axes class to which your figure should be added
    figsize: (width, height) or None | default = None
        optional parameter to define the figure size of the plot that is to be generated

    Returns
    -------
    None
        Figure is displayed
    """

    if adata.obs.get("n_counts") is None:
        adata.obs["n_counts"] = adata.X.sum(axis=1)

    if adata.obs.get("n_genes") is None:
        adata.obs["n_genes"] = np.sum(adata.X > 0, axis=1)

    x = adata.obs.get("n_counts").tolist()
    y = adata.obs.get("n_genes").tolist()

    data_plot = pd.DataFrame({"n_counts": x, "n_genes": y})
    min_UMI = 200

    ax = ax or plt.gca()

    ax.scatter(
        data=data_plot, x="n_counts", y="n_genes", alpha=0.4, s=0.5, rasterized=True
    )
    if figsize is not None:
        plt.figure(figsize=figsize)
    ax.hlines(max_genes, min_UMI, max(x), color="red", linestyle="dotted")
    ax.set_xlabel("n_counts")
    ax.set_ylabel("n_genes")
    ax.set_title("Filtering by the maximum gene count")

    return None


def max_mito(
    adata, max_mito=0.05, annotation_type="SYMBOL", species="human", copy=False, ax=None, figsize=None
):
    """visualize maximum mitochondrial gene percentage threshold.

    this function generates a knee-plot visualizing a given min_cells cutoff when given an adata object

    Parameters
    ----------
    adata: :class:`~anndata.AnnData`
        The annotated data matrix.
    threshold: `int` | default = 0
        integer value that defines the minimum expression threshold for a gene to be defined as expressed. Default value is 1.
    min_cells: `int` | default = 2
        visualize the chosen minimum number of cells that need to express a gene threshold
    ax: `axes` | default = None
        pass the axes class to which your figure should be added
    figsize: (width, height) or None | default = None
        optional parameter to define the figure size of the plot that is to be generated

    Returns
    -------
    None
        Figure is displayed

    """

    if (
        adata.obs.get("percent_mito") is not None
        and adata.obs.get("n_counts") is not None
    ):
        # generate data for plot
        x = adata.obs.get("n_counts").tolist()
        z = adata.obs.get("percent_mito").tolist()

        data_plot = pd.DataFrame({"n_counts": x, "percent_mito": z})

        # generate plot with cutoff
        ax = ax or plt.gca()

        ax.scatter(
            data=data_plot,
            x="n_counts",
            y="percent_mito",
            alpha=0.2,
            s=0.3,
            rasterized=True,
        )
        ax.set_xlabel("n_counts")
        ax.set_ylabel("percent_mito")
        ax.set_title("Mitochondrial gene content")
        if figsize is not None:
            plt.figure(figsize=figsize)
        ax.hlines(max_mito, 0, max(x), color="red", linestyle="dotted")

        return None

    elif annotation_type == "SYMBOL":
        if species == "human":

            # Extract mitochondrial genes
            mito_genes = [name for name in adata.var_names if name.startswith("MT-")]

            # for each cell compute fraction of counts in mito genes vs. all genes
            print("adding percent mitochondrial genes to dataframe for species human")
            adata.obs["percent_mito"] = (
                np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
            )

            # generate data for plot
            x = adata.obs.get("n_counts").tolist()
            z = adata.obs.get("percent_mito").tolist()

            data_plot = pd.DataFrame({"n_counts": x, "percent_mito": z})

            # generate plot with cutoff
            ax = ax or plt.gca()

            ax.scatter(
                data=data_plot,
                x="n_counts",
                y="percent_mito",
                alpha=0.2,
                s=0.3,
                rasterized=True,
            )
            ax.set_xlabel("n_counts")
            ax.set_ylabel("percent_mito")
            ax.set_title("mitochondrial gene content in dataset before filtering")
            if figsize is not None:
                plt.figure(figsize=figsize)
            ax.hlines(max_mito, 0, max(x), color="red", linestyle="dotted")

            if copy == True:
                return adata
            else:
                return None

        elif species == "mouse":
            # Extract mitochondrial genes
            mito_genes = [name for name in adata.var_names if name.startswith("mt-")]

            # for each cell compute fraction of counts in mito genes vs. all genes
            print("adding percent mitochondrial genes to dataframe for species mouse")
            adata.obs["percent_mito"] = (
                np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
            )

            # generate data for plot
            x = adata.obs.get("n_counts").tolist()
            z = adata.obs.get("percent_mito").tolist()

            data_plot = pd.DataFrame({"n_counts": x, "percent_mito": z})

            # generate plot with cutoff
            ax = ax or plt.gca()

            ax.scatter(
                data=data_plot,
                x="n_counts",
                y="percent_mito",
                alpha=0.2,
                s=0.3,
                rasterized=True,
            )
            ax.set_xlabel("n_counts")
            ax.set_ylabel("percent_mito")
            ax.set_title("mitochondrial gene content in dataset before filtering")
            if figsize is not None:
                plt.figure(figsize=figsize)
            ax.hlines(max_mito, 0, max(x), color="red", linestyle="dotted")

            if copy == True:
                return adata
            else:
                return None
        else:
            sys.exit(
                "currently only supports the species mouse and human. Please add implementation for additional species to besca."
            )

    elif annotation_type == "ENSEMBL":

        # read in reference mito file
        mito_list = get_mito_genes(species)
        mito_genes = [name for name in adata.var_names if name in mito_list]  # ensembl

        # for each cell compute fraction of counts in mito genes vs. all genes
        # the `.A1` is only necessary, as X is sparse - it transform to a dense array after summing
        adata.obs["percent_mito"] = (
            np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
        )
        # add the total counts per cell as observations-annotation to adata
        adata.obs["n_counts"] = adata.X.sum(axis=1).A1

        # generate data for plot
        x = adata.obs.get("n_counts").tolist()
        z = adata.obs.get("percent_mito").tolist()

        data_plot = pd.DataFrame({"n_counts": x, "percent_mito": z})

        # generate plot with cutoff
        ax = ax or plt.gca()
        ax.scatter(
            data=data_plot,
            x="n_counts",
            y="percent_mito",
            alpha=0.2,
            s=0.3,
            rasterized=True,
        )
        ax.set_xlabel("n_counts")
        ax.set_ylabel("percent_mito")
        ax.set_title("mitochondrial gene content in dataset before filtering")
        if figsize is not None:
            plt.figure(figsize=figsize)
        ax.hlines(max_mito, 0, max(x), color="red", linestyle="dotted")

        if copy == True:
            return adata
        else:
            return None

    else:
        print(
            'Calculating "percent_mito" and "n_genes" is currently only implemented for either SYMBOL gene symbols or ENSEMBL gene IDs.'
        )
        print(
            "Please supply an adata object where adata.obs contains these values already."
        )
