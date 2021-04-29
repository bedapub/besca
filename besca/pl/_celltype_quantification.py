import os
import sys
from matplotlib.pyplot import subplots, tight_layout
import seaborn as sns
import sys
from pandas import DataFrame, melt
from ..tl._count_occurrences import count_occurrence_subset, count_occurrence_subset_conditions

def celllabel_quant_boxplot(adata,
                            count_variable,
                            subset_variable,
                            condition_identifier,
                            plot_percentage = True,
                            condition_order=None,
                            myh=4,myw=8,mypal='Paired'):
    """ generate a box and whisker plot with overlayed swarm plot of celltype abundances

    This function takes a condition identifier to generate plots in which you can
    compare the celltype abundances between different conditions.

    Parameters
    ----------
    data: 'DataFrame'
        tidy dataframe that is outputed by the function besca.tl.count_occurrence_subset_conditions(..., make_tidy = True)
    count_variable: `str`
        obs on which to count - e.g. clusters or cell types
    subset_variable: `str`
        obs for stratification - e.g. donor, sample, experiment
    condition_identifier: `str`
        same condition_identifier used to generate data, used as hue - e.g. disease, storage, treatment
    plot_percentages: `bool` | default = True
        boolian indicating if the percentages or the total counts should be plotted
    condition_order: `list` | default = None
        list with the order in which the conditions should be plotted, for consistency
    myh: `integer` | default = 4
        height of the figure
    myw: `integer` | default = 8
        width of the figure
    mypal: `string` | default = "Paired"
        color palette for boxplots e.g. Paired, Blues_d        
    
    
    returns
    -------
    Figure

    """
    #generate dataframe to plot
    data = count_occurrence_subset_conditions(adata,
                                             subset_variable = subset_variable,
                                             count_variable = count_variable,
                                             condition_identifier = condition_identifier,
                                             return_percentage = plot_percentage)

    if condition_order==None:
        condition_order=list(adata.obs[condition_identifier].cat.categories)
        
    #add count_variable to dataframe as individual column
    data[count_variable] = data.index.tolist()

    #generate a tidy dataframe
    ylabel = 'percentage' if plot_percentage else 'count'
    data = melt(data, value_name = ylabel, id_vars = count_variable)

    #add columns to properly identify condition and subset
    new = data["variable"].str.split(" ", n = 2, expand = True)
    data[subset_variable] = new[1]
    data[condition_identifier] = new[2]
    data.drop(columns=['variable'], inplace=True)

    #create a new instance of a figure
    fig, (ax1) = subplots(1)
    fig.tight_layout()
    fig.set_figheight(myh)
    fig.set_figwidth(myw)

    #format the plot
    sns.set_style("white")
    sns.set_style("ticks")

    ax1 = sns.boxplot(data= data, x = count_variable, y = ylabel, hue = condition_identifier, palette=mypal, dodge = True, fliersize = 0,hue_order=condition_order)
    ax1 = sns.swarmplot(data= data, x = count_variable, y = ylabel, hue = condition_identifier, dodge=True, color='Black',hue_order=condition_order)
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles[:len(condition_order)], labels[:len(condition_order)])

    ax1.spines['bottom'].set_linewidth(1)
    ax1.spines['left'].set_linewidth(1)

    #shift the axis outwards and rotate labels
    sns.despine(offset=10, trim=True)
    ax1.tick_params(labelrotation=90,  length=6, width=2)

    #fix figure axis to include legend
    #tight_layout()

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
    >>> adata = bc.datasets.Kotliarov2020_processed()
    >>> adata.obs   = adata.obs.astype( {'batch' :  'category'})
    >>> fig = bc.pl.celllabel_quant_stackedbar(adata, count_variable = 'leiden', subset_variable = 'donor')

    .. plot::

        >>> import besca as bc
        >>> adata = bc.datasets.Kotliarov2020_processed()
        >>> adata.obs   = adata.obs.astype( {'batch' :  'category'})
        >>> fig = bc.pl.celllabel_quant_stackedbar(adata, count_variable = 'leiden', subset_variable = 'batch')
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

