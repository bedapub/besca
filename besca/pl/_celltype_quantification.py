import os
from matplotlib.pyplot import subplots, tight_layout
import seaborn as sns
from pandas import DataFrame, melt
from ..tl._count_occurrences import count_occurrence_subset, count_occurrence_subset_conditions

def celllabel_quant_boxplot(adata,
                            count_variable,
                            subset_variable,
                            condition_identifier,
                            plot_percentage = True):
    """ generate a box and whisker plot with overlayed swarm plot of celltype abundances

    This function takes a condition identifier to generate plots in which you can
    compare the celltype abundances between different conditions.

    Parameters
    ----------
    data: 'DataFrame'
        tidy dataframe that is outputed by the function besca.tl.compare_count_subsets(..., make_tidy = True)
    condition_identifier: `str`
        same condition_identifier used to generate data
    plot_percentages: `bool` | default = True
        boolian indicating if the percentages or the total counts should be plotted

    returns
    -------
    Figure

    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.pbmc_storage_processed_bbknn()
    >>> fig = bc.pl.celllabel_quant_boxplot(adata, count_variable = 'louvain', subset_variable = 'donor', condition_identifier = 'storage_condition',  plot_percentage = True);

    .. plot::

        >>> import besca as bc
        >>> adata = bc.datasets.pbmc_storage_processed_bbknn()
        >>> fig = bc.pl.celllabel_quant_boxplot(adata, count_variable = 'louvain', subset_variable = 'donor', condition_identifier = 'storage_condition',  plot_percentage = True);

    """
    #generate dataframe to plot
    if plot_percentage:
        data = count_occurrence_subset_conditions(adata, 
                                                 subset_variable = subset_variable, 
                                                 count_variable = count_variable, 
                                                 condition_identifier = condition_identifier,  
                                                 return_percentage = True)
    else:
        data = count_occurrence_subset_conditions(adata, 
                                                 subset_variable = subset_variable, 
                                                 count_variable = count_variable, 
                                                 condition_identifier = condition_identifier,  
                                                 return_percentage = False)

    #add count_variable to dataframe as individual column
    data[count_variable] = data.index.tolist()

    #generate a tidy dataframe
    if plot_percentage:
        data = melt(data, value_name = 'percentage', id_vars = count_variable)
    else:
        data = melt(data, value_name = 'count', id_vars = count_variable)

    #add columns to properly identify condition and subset
    new = data["variable"].str.split(" ", n = 2, expand = True)
    data[subset_variable] = new[1]
    data[condition_identifier] = new[2]
    data.drop(columns=['variable'], inplace=True)

    #create a new instance of a figure
    fig, (ax1) = subplots(1)
    fig.tight_layout()
    fig.set_figheight(4)
    fig.set_figwidth(8)

    #format the plot
    sns.set_style("white")
    sns.set_style("ticks")

    if plot_percentage:
        ax1 = sns.boxplot(data= data, x = count_variable, y = 'percentage', hue = condition_identifier, palette="Blues_d", dodge = True, fliersize = 0)
        ax1 = sns.swarmplot(data= data, x = count_variable, y = 'percentage', hue = condition_identifier, dodge=True, color='Black')
    else: 
        ax1 = sns.boxplot(data= data, x = count_variable, y = 'count', hue = condition_identifier, palette="Blues_d", dodge = True, fliersize = 0)
        ax1 = sns.swarmplot(data= data, x = count_variable, y = 'count', hue = condition_identifier, dodge=True, color='Black')

    ax1.spines['bottom'].set_linewidth(1)
    ax1.spines['left'].set_linewidth(1)

    #shift the axis outwards and rotate labels
    sns.despine(offset=10, trim=True)
    ax1.tick_params(labelrotation=90,  length=6, width=2)

    #fix figure axis to include legend
    tight_layout()

    #return the generated figure
    return(fig)

def celllabel_quant_stackedbar(adata,
                               subset_variable,
                               count_variable = 'celltype'):
    """ Generate a stacked bar plot of the percentage of labelcounts within each AnnData subset

    parameters
    ----------
    adata: AnnData
      the AnnData object 
    subset_variable: `str`
        string identifying the column name of adata.obs along which the AnnData object should be subsetted
    count_variable: `str` | default  = 'celltype'
        string identiyfing the column of adata.obs containing the labels to be counted
    
    returns
    -------
    Figure
    
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.pbmc_storage_processed_bbknn()
    >>> fig = bc.pl.celllabel_quant_stackedbar(adata, count_variable = 'louvain', subset_variable = 'donor');

    .. plot::

        >>> import besca as bc
        >>> adata = bc.datasets.pbmc_storage_processed_bbknn()
        >>> fig = bc.pl.celllabel_quant_stackedbar(adata, count_variable = 'louvain', subset_variable = 'donor');
    
    """
    counts_celltype = count_occurrence_subset(adata, subset_variable, count_variable =count_variable, return_percentage = False)
    percentages = DataFrame(index = counts_celltype.index.tolist(), columns = counts_celltype.columns.tolist())
    
    for cell in counts_celltype.index.tolist():
        data = counts_celltype.loc[cell].tolist()
        percentage = [x/sum(data) for x in data]
        percentages.loc[cell] = percentage

    fig = percentages.plot(kind='bar', stacked=True, figsize=(8, 4))
    fig.set_ylabel('percentage')
    fig.legend(bbox_to_anchor=(1, 1))
    
    #fix figure axis to include legend
    tight_layout()

    #return figure
    return(fig)

def louvain_quant_stackedbar(adata,
                             subset_variable,
                             count_variable = 'louvain'):
    """ Generate a stacked bar plot of the percentage of cells in each louvain cluster from each AnnData subset

    In contrast to the celllabel_quant_stackedbar function which can also be used to visualize the lovain clusters
    this function also sorts the louvain clusters in proper order on the x-Axis.

    parameters
    ----------
    adata: AnnData
      the AnnData object 
    subset_variable: `str`
        string identifying the column name of adata.obs along which the AnnData object should be subsetted
    count_variable: `str` | default  = 'celltype'
        string identiyfing the column of adata.obs containing the labels to be counted
    
    returns
    -------
    Figure

    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.pbmc_storage_processed_bbknn()
    >>> fig = bc.pl.louvain_quant_stackedbar(adata, subset_variable = 'donor');

    .. plot::

        >>> import besca as bc
        >>> adata = bc.datasets.pbmc_storage_processed_bbknn()
        >>> fig = bc.pl.louvain_quant_stackedbar(adata, subset_variable = 'donor');
    
    """
    counts_celltype = count_occurrence_subset(adata, subset_variable, count_variable =count_variable, return_percentage = False)
    percentages = DataFrame(index = counts_celltype.index.tolist(), columns = counts_celltype.columns.tolist())
    
    for cell in counts_celltype.index.tolist():
        data = counts_celltype.loc[cell].tolist()
        percentage = [x/sum(data) for x in data]
        percentages.loc[cell] = percentage

    percentages.index = list(map(int, percentages.index.tolist()))
    percentages = percentages.sort_index()

    fig = percentages.plot(kind='bar', stacked=True, figsize=(8, 4))
    fig.set_ylabel('percentage')
    fig.set_xlabel('louvain cluster')
    fig.legend(bbox_to_anchor=(1, 1))

    #fix figure axis to include legend
    tight_layout()

    #return figure
    return(fig)