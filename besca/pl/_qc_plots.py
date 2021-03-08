# Quality control plots

from pandas import DataFrame
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy
from sklearn import linear_model
from numpy import ndarray

def transcript_capture_efficiency(adata, ax = None):
    """Plot total gene counts vs detection probability.

    Visualize the transcript capture efficiency curve
    by plotting the total gene counts on the x-axis and
    the detection probability on the y-axis.

    parameters
    ----------
    adata: `AnnData`
        AnnData object containing data that is to be visualized. Genes need to be in coloumns and cells in rows
    ax: `axes` | default = None
        pass the axes class to which your figure should be added
    
    returns
    -------
    None 
        the figure is displayed

    Example
    -------
    Display transcript capture efficiency plot.
    
    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> adata = bc.datasets.pbmc3k_raw()
    >>> fig, ax = plt.subplots(1)
    >>> bc.pl.transcript_capture_efficiency(adata,ax=ax)
    
    .. plot:: 

        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> adata = bc.datasets.pbmc3k_raw()
        >>> fig, ax = plt.subplots(1)
        >>> bc.pl.transcript_capture_efficiency(adata,ax=ax)

    """
    #get adata object

    if type(adata.X) == np.ndarray:
        adata_X = adata.X
    else:
        adata_X = adata.X.toarray()

    #extract information to be plotted
    fraction_pos = sum(adata_X > 0) / adata.shape[0]
    total_gene_counts = np.log2(sum(adata_X)+1)

    #generate dataframe containing information
    data = DataFrame(data={'fraction_pos':fraction_pos, 'total_gene_counts':total_gene_counts})

    #generate figure
    ax = ax or plt.gca()

    ax.scatter(data=data, x='total_gene_counts', y = 'fraction_pos', alpha = 0.4, s = 4, edgecolors='none', rasterized = True)
    ax.set_title('transcript capture efficiency')
    ax.set_xlabel('log2(total gene counts)')
    ax.set_ylabel('detection probability')

def library_size(adata,
                 ax = None,
                 bins = 100):
    """ Plot library size.

    Generates a histogram of the library size per cell.

    parameters
    ----------
    adata: `AnnData`
        AnnData object containing data that is to be visualized. Genes need to be in coloumns and cells in rows
    ax: `axes` | default = None
        pass the axes class to which your figure should be added
    bins: `int` | default = 100
        the number of bins that should be shown on the histrogram
    
    returns
    -------
    None 
        displays a figure

    Example
    -------
    Plot distribution of librarysize from an example dataset.
    
    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> adata = bc.datasets.pbmc3k_raw()
    >>> fig, ax = plt.subplots(1)
    >>> bc.pl.library_size(adata,ax=ax)
    
    .. plot:: 

        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> adata = bc.datasets.pbmc3k_raw()
        >>> fig, ax = plt.subplots(1)
        >>> bc.pl.library_size(adata,ax=ax)

    """
    if type(adata.X) == np.ndarray:
        adata_X = adata.X
    else:
        adata_X = adata.X.toarray()

    library_size = np.sum(adata_X, axis = 1)/1e6

    ax = ax or plt.gca()

    sns.distplot(library_size, ax = ax, bins = bins);
    ax.set_title('library size distribution')
    ax.set_ylabel('number of cells')
    ax.set_xlabel('library size (millions)')

def detected_genes(adata,
                   ax = None,
                   bins = 100):
    """ Plot number of detected genes.

    Generates a histogram of the number of detected genes per cell.
    
    parameters
    ----------
        adata: `AnnData`
        AnnData object containing data that is to be visualized. Genes need to be in coloumns and cells in rows
    ax: `axes` | default = None
        pass the axes class to which your figure should be added
    bins: `int` | default = 100
        the number of bins that should be shown on the histrogram
    
    returns
    -------
    None 
        displays the figure

    Example
    -------
    Plot number of detected genes from an example dataset.
    
    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> adata = bc.datasets.pbmc3k_raw()
    >>> fig, ax = plt.subplots(1)
    >>> bc.pl.detected_genes(adata,ax=ax)
    
    .. plot:: 

        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> adata = bc.datasets.pbmc3k_raw()
        >>> fig, ax = plt.subplots(1)
        >>> bc.pl.detected_genes(adata,ax=ax)

    """
    if type(adata.X) == np.ndarray:
        adata_X = adata.X
    else:
        adata_X = adata.X.toarray()

    NODG = np.sum(adata_X > 0, axis = 1)
    
    sns.distplot(NODG, ax = ax, kde = False, norm_hist = False, bins = bins);
    ax.set_title('NODG');
    ax.set_ylabel('number of cells');
    ax.set_xlabel('number of detected genes');

def dropouts(adata,
             ax = None,
             bins = 100):
    """ Plot number of dropouts.

    Generates a histrogram showing the number of dropouts per cell.
    
    parameters
    ----------
    adata: `AnnData`
        AnnData object containing data that is to be visualized. Genes need to be in coloumns and cells in rows
    ax: `axes` | default = None
        pass the axes class to which your figure should be added
    bins: `int` | default = 100
        the number of bins that should be shown on the histrogram
    
    returns
    -------
    None 
        displays the figure

    Example
    -------
    Plot number of dropout genes from an example dataset.
    
    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> adata = bc.datasets.pbmc3k_raw()
    >>> fig, ax = plt.subplots(1)
    >>> bc.pl.dropouts(adata,ax=ax)
    
    .. plot:: 

        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> adata = bc.datasets.pbmc3k_raw()
        >>> fig, ax = plt.subplots(1)
        >>> bc.pl.dropouts(adata,ax=ax)   

    """
    if type(adata.X) == np.ndarray:
        adata_X = adata.X
    else:
        adata_X = adata.X.toarray()

    dropouts = np.sum(adata_X == 0, axis = 1)

    sns.distplot(dropouts, ax=ax, kde=False, norm_hist=False, bins=bins)
    ax.set_title('dropouts')
    ax.set_xlabel('number of dropouts')
    ax.set_ylabel('number of cells')

def librarysize_overview(adata,
                         bins=100):
    """ Generates overview figure of libarysize, dropouts and detected genes.

    Generates one overview figure showing histograms of the librarysize, number of
    detected genes and number of dropouts per cell. Can be used instead of plotting
    librarysize(), dropouts() abd detected_genes() seperately.
    
    parameters
    ----------
    adata: `AnnData`
        AnnData object containing data that is to be visualized. Genes need to be in coloumns and cells in rows
    bins: `int` | default = 100
        the number of bins that should be shown on the histrogram
    
    returns
    -------
    Figure 
        returns a figure (is also displayed)

    Example
    -------

    Generate overview of the characterisitcs of an example dataset.
    
    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> adata = bc.datasets.pbmc3k_raw()
    >>> bc.pl.librarysize_overview(adata)
    
    .. plot:: 

        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> adata = bc.datasets.pbmc3k_raw()
        >>> bc.pl.librarysize_overview(adata)
        
    """
    if type(adata.X) == np.ndarray:
        adata_X = adata.X
    else:
        adata_X = adata.X.toarray()

    library_size = np.sum(adata_X, axis = 1)/1e6
    dropouts = np.sum(adata_X == 0, axis = 1)
    NODG = np.sum(adata_X > 0, axis = 1)

    fig = plt.figure(figsize=(11, 8)) 
    ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=2)
    ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=1)
    ax3 = plt.subplot2grid((2, 2), (1, 1), colspan=1, sharey=ax2)

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.1, hspace = 0.3)

    #plot librarysize
    sns.distplot(library_size, ax = ax1, bins = bins);
    ax1.set_title('library size distribution')
    ax1.set_ylabel('number of cells')
    ax1.set_xlabel('library size (millions)')
    
    sns.set_style("white")
    sns.set_style("ticks")
    
    sns.distplot(NODG, ax = ax2, kde = False, norm_hist = False, bins = bins);
    ax2.set_title('NODG');
    ax2.set_ylabel('number of cells');
    ax2.set_xlabel('number of detected genes', labelpad=16);
    
    sns.despine(offset=10, trim=True, ax = ax2)
    
    sns.distplot(dropouts, ax = ax3, kde= False, norm_hist = False, bins = bins);
    ax3.set_title('dropouts');
    ax3.set_xlabel('number of dropouts', labelpad = 16);
    
    sns.despine(offset=10, trim=True, ax = ax3)
    return(fig)


def top_genes_counts(adata, top_n=25, ax=None):
    """ plot top n genes that contribute to fraction of counts per cell

    Generate box and whisker plot of the fraction of counts in a cell per gene
    which shows the top n genes that most strongly contribute to UMI counts.

    parameters
    ----------
    adata: `AnnData`
        AnnData object containing data that is to be visualized.
        Genes need to be in columns and cells in rows
    top_n: `int` | default = 25
        number of genes that should be visualized in the plot
    ax: `axes` | default = None
        the axes instance to which your plot should be added

    returns
    -------
    None
        If an ax object is passed to the function the function returns nothing
    Figure
        If ax = None then a figure is returned

    Example
    -------

    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> adata = bc.datasets.pbmc3k_raw()
    >>> bc.pl.top_genes_counts(adata)

    .. plot::

        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> adata = bc.datasets.pbmc3k_raw()
        >>> bc.pl.top_genes_counts(adata)

    """
    # calculate total counts and frac_reads
    adata.var['total_counts'] = sum(adata.X.toarray())
    adata.var['frac_reads'] = sum(adata.X.toarray()) / sum(sum(adata.X.toarray()))
    if type(adata.X) == ndarray:
        adata_X = adata.X
    else:
        adata_X = adata.X.toarray()

    # calculate counts per cell
    counts_per_cell = adata_X.sum(axis=1)

    # calculate fraction of counts originating from each gene  per cell
    fracs = adata_X / counts_per_cell[:, None]

    # make dataframe to sort
    data = DataFrame(data=fracs, columns=adata.var_names).T
    data['total_counts'] = adata.var.total_counts
    data['frac_reads'] = adata.var.frac_reads
    data.sort_values(ascending=False, by='total_counts', inplace=True)

    # get top n values
    data = data.head(top_n)

    # calculate cumsum
    cum_sum = sum(data.frac_reads)

    # remove unwanted columns for plotting
    data.drop(columns=['total_counts', 'frac_reads'], inplace=True)

    # plot results
    data = data.T

    # define the properties of the outliermarkers so they are not as prominent
    flierprops = dict(marker=".", markersize=2, linestyle='none')

    if ax is not None:
        ax = sns.boxplot(data=data, orient="h", color='orange', width=0.7,
                         flierprops=flierprops)
        ax.set_ylabel('gene')
        ax.set_xlabel('fraction of UMI counts per cell')
        ax.set_title('Top ' + str(top_n) + ' genes account for ' +
                     str(round(cum_sum*100, 2)) + '% of all UMI counts')
        return(None)
    else:
        # generate new figure instance that is returned
        fig = plt.figure(figsize=(12/25*top_n, 6))
        fig = sns.boxplot(data=data, orient="h", color='orange', width=0.7,
                          flierprops=flierprops)
        fig.set_ylabel('gene')
        fig.set_xlabel('fraction of UMI counts per cell')
        fig.set_title('Top ' + str(top_n) + ' genes account for ' +
                      str(round(cum_sum*100, 2)) + '% of all UMI counts')
        return(fig)
