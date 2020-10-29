import seaborn as sns
from matplotlib.pyplot import figure, subplots_adjust, tight_layout
import sys

def split_violin(tidy_data,
                 x_axis,
                 y_axis,
                 split_variable,
                 order = None,
                 ax = None,
                 inner = 'box'):
    """plot ssplit violin plots.
    
    General plotting function to produce split violin plots.

    parameters
    ----------
    tidy_data: DataFrame
        pandas DataFrame containing the complete data that is to be plotted in a tidy format.
    x_axis: `str`
        string identifying which column of the DataFrame is to be plotted on the x-axis
    y_axis: `str`
        string identifying which column of the DataFrame is to be plotted on the y-axis
    split_variable: `str`
        string identifying which column of the DataFrame is to be used to generate the 
        split violin plot (can only contain two categories of data!)
    subset_variable_label: `str`
        string identifiyng which column of the DataFrame contains the variables that 
        should be used to make datasubsets for each plot of the stacked violin plot
    subset_variabels: `list`
        list identifying the subsets that should be generated
    fig_width: `int` | default = 8
        int value defining figure width
    fig_height: `int` | default = 4
        int value defining figure height of one figure
    ax: `axes` | default = None
        pass the axes class to which your figure should be added, if none is supplied a new figure is generated

    returns
    -------
    Figure
    """
    ax = ax or plt.gca()

    #set plotting style
    sns.set_style("white")
    sns.set_style("ticks")

    ax = sns.violinplot(x=x_axis, y=y_axis, hue=split_variable,  data=tidy_data, palette="muted", split=True, order = order, inner = inner)

    #shift spines outwards
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)

    #remove spines
    sns.despine(offset=10, trim=True)

    #make labels 90 degrees so they are readible
    ax.tick_params(labelrotation=90,  length=6, width=2)

    return(None)

def box_per_ind(plotdata,y_axis,x_axis,order=None,fig_width = 4,fig_height = 3):
    """plot boxplot with values per individual.
    
    General plotting function to produce one or multiple boxplots for 
    average/fraction gene expression per individual/sample.

    parameters
    ----------
    plotdata: DataFrame
        pandas DataFrame containing the complete data that is to be plotted.
    x_axis: `str`
        string identifying which column of the DataFrame is to be plotted on the x-axis (condition)
    y_axis: `list` or `str`
        string identifying which column of the DataFrame is to be plotted on the y-axis (genes)
    order: `list`
        list identifying the order for the categories on the x_axis
    fig_width: `int` | default = 4
        int value defining figure width
    fig_height: `int` | default = 3
        int value defining figure height of one figure
        
    returns
    -------
    Figure
    """
    
    #set plotting style
    sns.set_style("white")
    sns.set_style("ticks")

    
    if isinstance(y_axis, list)==False:
        y_axis=[y_axis]
    
    #only retain elements present
    y_axisk=[]
    for y in y_axis:
        y_axisk=y_axisk+list(set(y).intersection(set(plotdata.columns)))
    y_axis=y_axisk
        
    #check that we have genes
    if len(y_axis)==0:
        sys.exit('Please select valid gene names')

    #check that we have the condition
    if (x_axis in plotdata.columns)==False:
        sys.exit('Please select a valid condition name')
        
    if order==None:
        order=list(set(plotdata[x_axis]))
        
    #determine number of subplots
    number_of_subplots=len(y_axis)
    
    #initiate figure
    fig = figure()
    
    #adjust size of figure if desired
    if fig_width is not None:
        fig.set_figwidth(fig_width)
    if fig_height is not None:
        fig.set_figheight(fig_height * number_of_subplots)
    
    #adjust amount of space between subplots
    subplots_adjust(hspace=0.000)

    #########################################################
    #plot first figure (this adds the legend above the plot)
    #########################################################
    
    #plot figure
    ax0 = fig.add_subplot(number_of_subplots,1,1)
    ax0 = sns.boxplot(x=x_axis, y=y_axis[0], data=plotdata, order=order,palette="muted")
    ax0 = sns.stripplot(x=x_axis, y=y_axis[0], data=plotdata, order=order,color='black')
    ax0.axes.get_xaxis().set_visible(False)
    ax0.yaxis.tick_right()

    #get correct label for the y-axis
    ax0.set_ylabel(y_axis[0])

    #move legend above the plot
    #ax0.legend(loc=9, bbox_to_anchor=(0.5, 1.5), ncol=2)


    #########################################################
    #plot all subsequent figures dynamically
    #########################################################    
    if number_of_subplots >=2:
        
        for i in range(1, number_of_subplots):
            
            #get indicator for subplot number
            v = i+1
            
            #add subplot
            ax1 = fig.add_subplot(number_of_subplots,1,v, sharey = ax0)
            ax1 = sns.boxplot(x=x_axis, y=y_axis[i], data=plotdata, order=order,palette="muted")
            ax1 = sns.stripplot(x=x_axis, y=y_axis[i], data=plotdata, order=order,color='black')
            ax1.axes.get_xaxis().set_visible(False)
            ax1.yaxis.tick_right()

            #get correct label for the y-axis
            ax1.set_ylabel(y_axis[i])

            #remove legend since we only need it once
            #ax1.get_legend().remove()

        #set x-axis on the last plot generated to visible
        ax1.axes.get_xaxis().set_visible(True)
        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)
    
    else:
        ax0.axes.get_xaxis().set_visible(True)
        ax0.set_xticklabels(ax0.get_xticklabels(), rotation=90)

    tight_layout()
    subplots_adjust(hspace=0.000)


    return(None)

def stacked_split_violin(tidy_data,
                         x_axis,
                         y_axis,
                         split_variable,
                         subset_variable_label,
                         subset_variables,
                         fig_width = 8,
                         fig_height = 4,
                         order = None,
                         inner = 'box'):
    """ plot stacked split violin plots.
    
    General plotting function to produce stacked split violin plots.

    parameters
    ----------
    tidy_data: DataFrame
        pandas DataFrame containing the complete data that is to be plotted in a tidy format.
    x_axis: `str`
        string identifying which column of the DataFrame is to be plotted on the x-axis
    y_axis: `str`
        string identifying which column of the DataFrame is to be plotted on the y-axis
    split_variable: `str`
        string identifying which column of the DataFrame is to be used to generate the 
        split violin plot (can only contain two categories of data!)
    subset_variable_label: `str`
        string identifiyng which column of the DataFrame contains the variables that 
        should be used to make datasubsets for each plot of the stacked violin plot
    subset_variabels: `list`
        list identifying the subsets that should be generated
    fig_width: `int` | default = 8
        int value defining figure width
    fig_height: `int` | default = 4
        int value defining figure height of one figure
    order:
    inner: 'box' or 'quartile' or 'point' or 'stick'
        define how the datapoints should be displayed in the violin interior, see seaborns documentation for more details

    returns
    -------
    Figure
        
    """
    
    #determine number of subplots
    number_of_subplots=len(subset_variables)
    
    #initiate figure
    fig = figure()
    
    #adjust size of figure if desired
    if fig_width is not None:
        fig.set_figwidth(fig_width)
    if fig_height is not None:
        fig.set_figheight(fig_height * number_of_subplots)
    
    #adjust amount of space between subplots
    subplots_adjust(hspace=0.000)
    
    #########################################################
    #plot first figure (this adds the legend above the plot)
    #########################################################
    
    #generate datasubset
    data = tidy_data[tidy_data.get(subset_variable_label) == subset_variables[0]];
    
    #plot figure
    ax0 = fig.add_subplot(number_of_subplots,1,1)
    ax0 = sns.violinplot(x=x_axis, y=y_axis, hue=split_variable,  data=data, palette="muted", split=True, order = order, inner = inner)
    ax0.axes.get_xaxis().set_visible(False)
    ax0.yaxis.tick_right()

    
    #get correct label for the y-axis
    ax0.set_ylabel(subset_variables[0])

    #move legend above the plot
    ax0.legend(loc=9, bbox_to_anchor=(0.5, 1.5), ncol=2)
    
    #########################################################
    #plot all subsequent figures dynamically
    #########################################################    
    if number_of_subplots >=2:
        
        for i in range(1, number_of_subplots):
            #subset data
            data = tidy_data[tidy_data.get(subset_variable_label) == subset_variables[i]]
            
            #get indicator for subplot number
            v = i+1
            
            #add subplot
            ax1 = fig.add_subplot(number_of_subplots,1,v, sharey = ax0)
            ax1 = sns.violinplot(x=x_axis, y=y_axis, hue=split_variable,  data=data, palette="muted", split=True, order = order, inner = inner)
            ax1.axes.get_xaxis().set_visible(False)
            ax1.yaxis.tick_right()

            #get correct label for the y-axis
            ax1.set_ylabel(subset_variables[i])

            #remove legend since we only need it once
            ax1.get_legend().remove()

        #set x-axis on the last plot generated to visible
        ax1.axes.get_xaxis().set_visible(True)
        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)
    
    else:
        ax0.axes.get_xaxis().set_visible(True)
        ax0.set_xticklabels(ax0.get_xticklabels(), rotation=90)

    tight_layout()
    subplots_adjust(hspace=0.000)
    return(fig)