# generate split plots too compare gene expression between two cellpopulations
import seaborn as sns
from pandas import concat, melt, DataFrame
from matplotlib.pyplot import subplots, gca
from matplotlib.offsetbox import AnchoredText
import sys
import os
from besca.pl._general import stacked_split_violin, split_violin
from besca._helper import get_raw
import numpy as np


def _make_tidy(adata, genes, split_variable, label_split_variable=None, use_raw=True):
    """Helper function to generate the tidy dataframes necessary to plot split_gene expression.

    This helper function should only be called in the context of the split gene expression plots.

    parameters
    ----------
    adata: `AnnData`
        AnnData object that contains the data you wish to plot
    genes: `[str]`
        list of strings identifiyng the genes which are to be plotted
    split_variable: `str`
        string identfying the column which contains the annotation which cell belongs to which group (should only
        contain 2 different groups!!) If you have more than 2 groups you will need to look into other visualization techniques.
    label_split_variable: `str` | default = None
        string indicating which which label the legend should be displayed
    use_raw: `bool` | default = True
        boolian variable indicating if the gene expression values saved in adata.raw should be plotted or not
    """

    # get labels if not previously defined
    if label_split_variable is None:
        label_split_variable = split_variable
    # get raw object if using
    bdata = get_raw(adata) if use_raw else adata.copy()

    # convert X to array
    if type(bdata.X) == np.matrixlib.defmatrix.matrix:
        pass
    else:
        bdata.X = bdata.X.todense()

    # get index of raw gene names
    gene_index = bdata.var_names.tolist()

    # get index location of genes
    index = []
    for gene in genes:
        index.append(gene_index.index(gene))

    # extract AnnData object containing only the relevant genes (this decreases computational time later)
    bdata = bdata[:, index]
    # check that group variable only contains two groups
    groups = bdata.obs.get(split_variable).value_counts().index.tolist()
    # perform some basic sanity checks
    if len(groups) > 2:
        sys.exit(
            "please ensure that the groups you wish to plot againt one another are not more than 2!"
        )
    if len(groups) == 1:
        group1 = groups[0]
        # extract matrix of gene expression values for each group
        matrix_group1 = bdata[bdata.obs.get(split_variable) == group1].X
        matrix_group2 = None
    elif len(groups) == 2:
        # get names of two groups we are plotting against each other
        group1 = groups[0]
        group2 = groups[1]
        # extract matrix of gene expression values for each group
        matrix_group1 = bdata[bdata.obs.get(split_variable) == group1].X
        matrix_group2 = bdata[bdata.obs.get(split_variable) == group2].X
    else:
        print("no groups")
    # check that the number of cells contained in the subset is larger than 1
    if len(matrix_group1.shape) >= 1:
        condition1 = True
    else:
        condition1 = None
    if matrix_group2 is not None:
        if len(matrix_group2.shape) >= 1:
            condition2 = True
        else:
            condition2 = None
    else:
        condition2 = None
    # generate tidy frame for datasubset of group1
    if condition1:
        # initialize an empty dataframe to save our results for each group into
        gene_expression_group1 = DataFrame(index=range(matrix_group1.shape[0]))
        if len(matrix_group1.shape) > 1:
            # extract expression value for each of our genes
            # print('iterating through genes for group1')
            for i in range(0, len(genes)):
                # get values for group1
                gene_expression_group1[genes[i]] = matrix_group1[:, i]
        else:
            # extract expression value for each of our genes
            # print('iterating through genes for group1')
            for i in range(0, len(genes)):
                # get values for group1
                gene_expression_group1[genes[i]] = matrix_group1

        # make into tidy dataframe and add annotation
        gene_expression_group1 = melt(
            gene_expression_group1, value_name="expression", var_name="gene"
        )
        gene_expression_group1[label_split_variable] = group1

    else:
        print("The dataset containts no cells for condition ", group2)
    # generate tidy frame for datasubset of group2
    if condition2:
        # initialize an empty dataframe to save our results for each group into
        gene_expression_group2 = DataFrame(index=range(matrix_group2.shape[0]))
        if len(matrix_group2.shape) > 1:
            # extract expression value for each of our genes
            # print('iterating through genes for group1')
            for i in range(0, len(genes)):
                # get values for group1
                gene_expression_group2[genes[i]] = matrix_group2[:, i]
        else:
            # extract expression value for each of our genes
            # print('iterating through genes for group1')
            for i in range(0, len(genes)):
                # get values for group1
                gene_expression_group2[genes[i]] = matrix_group2
        # make into tidy dataframe and add annotation
        gene_expression_group2 = melt(
            gene_expression_group2, value_name="expression", var_name="gene"
        )
        gene_expression_group2[label_split_variable] = group2
    else:
        print("The dataset contains no cells for second condition.")
    # merge dataframes from the two conditions into one
    if condition1 is not None and condition2 is not None:
        # merge two pd dataframes
        data = concat([gene_expression_group1, gene_expression_group2], axis=0)
    elif condition1 is not None:
        data = gene_expression_group1
    elif condition2 is not None:
        data = gene_expression_group2
    # return data and successfully exit function
    return data
    sys.exit(0)


def gene_expr_split(
    adata, genes, split_variable=None, label_split_variable=None, use_raw=True, ax=None
):
    """visualize gene expression of two groups as a split violin plot

    This function generates split violin plots to compare the gene expression levels of
    a subset of genes between two groups.

    If a split_variable is provided it first subsets the provided adata object according to this variable
    and generates an individual plot for each of the subsets. In the default run configuration, this function
    generates split violin plots where the x-axis contains the genes of interest, the y-axis the expression values
    and the hue of the violin plot represents the two groups that are being compared. If plot_genes = True, then
    this function returns an individual plot for each gene. Here if a split_variable is used to subset the data the
    x-axis then contains all the different datasubsets.

    If using this function to subset data based on the split_variable the generated figures must be written to an
    output directory. The same applies for generating one plot per gene.

    If only looking at one data subset, the figure can also be returned as an object if outdir = None.

    parameters
    ----------
    adata : :class:`~anndata.AnnData`
        The annotated data matrix used as input
    genes: `list`
        List of gene ids for which the gene expression levels should be compared
    split_variable: `str`
        string identfying the column which contains the annotation which cell belongs to which group (should only
        contain 2 different groups!!) If you have more than 2 groups you will need to look into other visualization techniques.
    use_raw: `bool` (default = True)
        boolian indicator if adata.raw data should be used as input for gene expression levels
    plot_genes: `bool` | default = False
        boolian indicator if one individual plot should be generated per gene or if all the genes should be on the same plot
    outdir: `str`
        filepath indicating where the generated figures should be saved
    fig_width: `int` | default = 15
        width of the generated figure
    fig_height: `int` | default = 8
        height of the generated figure
    dpi: `int`
        dpi resolution of figures that are saved to file
    name: `str` (default = 'compare_gene_expression')
        name of generated plots and prefix to exported files

    returns
    -------
    None
        writes figures to file specified in outdir.

    """
    # get labels if not previously defined
    if label_split_variable is None:
        label_split_variable = split_variable
    ### get tidy dataframe
    data_merged = _make_tidy(
        adata=adata,
        genes=genes,
        split_variable=split_variable,
        label_split_variable=label_split_variable,
        use_raw=use_raw,
    )

    ax = ax or gca()

    split_violin(
        tidy_data=data_merged,
        x_axis="gene",
        y_axis="expression",
        split_variable=label_split_variable,
        order=None,
        ax=ax,
    )
    return None


def gene_expr_split_stacked(
    adata,
    genes,
    split_variable,
    label_split_variable=None,
    subset_variable="celltype",
    label_subset_variable=None,
    use_raw=True,
    fig_width=None,
    fig_height=None,
    order=None,
    inner="stick",
):

    """
    Stacked violin plot for visualization of genes expression.

    parameters
    ----------
    adata : :class:`~anndata.AnnData`
        The annotated data matrix used as input
    genes: `list`
        List of gene ids for which the gene expression levels should be compared
    split_variable: `str`
        string identfying the column which contains the annotation which cell belongs to which group (should only
        contain 2 different groups!!) If you have more than 2 groups you will need to look into other visualization techniques.
    subset_variable: `str`
        string identifying the column along which the AnnData object should be split into subsets. Can contain as many groups as you wish.
        parameters
    use_raw: `bool` (default = True)
        boolian indicator if adata.raw data should be used as input for gene expression levels
    fig_width: `int` | default = 15
        width of the generated figure
    fig_height: `int` | default = 8
        height of the generated figure
    order: lists of strings
        Order to plot the categorical levels in
    inner: `str` (default = 'stick')
        see seaborn violin plot.

    returns
    -------
    fig

    """
    if label_subset_variable is not None:
        label_subset_variable = label_subset_variable
    else:
        label_subset_variable = subset_variable

    if label_split_variable is not None:
        label_split_variable = label_split_variable
    else:
        label_split_variable = split_variable

    # get tidy data

    # define subsets
    subsets = adata.obs.get(subset_variable).value_counts().index.tolist()

    data_to_merge = {}
    for subset in subsets:
        adata_subset = adata[adata.obs.get(subset_variable) == subset]

        data = _make_tidy(
            adata=adata_subset,
            genes=genes,
            split_variable=split_variable,
            label_split_variable=label_split_variable,
            use_raw=use_raw,
        )
        data[label_subset_variable] = subset

        data_to_merge.update({subset: data})

    # merge all iterations into one dataframe
    print("merging a total of ", str(len(subsets)), " datasubset")

    # print('merging datasubset 1')
    data_merged = data_to_merge.get(subsets[0])
    for iteration in subsets[1:]:
        # print('merging subset' + str(variables.index(iteration) + 1))
        data_merged = concat([data_merged, data_to_merge.get(iteration)], axis=0)

    fig = stacked_split_violin(
        tidy_data=data_merged,
        x_axis=label_subset_variable,
        y_axis="expression",
        split_variable=label_split_variable,
        subset_variable_label="gene",
        subset_variables=genes,
        fig_width=fig_width,
        fig_height=fig_height,
        order=order,
        inner=inner,
    )
    return fig
