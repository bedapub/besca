import sys
from matplotlib.patches import Circle
from matplotlib.colors import Normalize
from math import log10, ceil
from matplotlib.pyplot import subplots, get_cmap, figure
from matplotlib import gridspec
from matplotlib.colorbar import ColorbarBase
from pandas import DataFrame, Series
from numpy import ndarray, arange, float32

# import helper functions from besca
from besca._helper import get_raw

import pytest


def _generate_circle(expression_values, center, radius, ax):
    """helper function to plot heatmap circle with black outline"""
    slices, labels = ax.pie(
        expression_values.counts,
        colors=expression_values.color,
        shadow=False,
        center=center,
        radius=radius,
        startangle=90,
    )

    # properly set wedge edge color
    for i in expression_values.index:
        slices[i].set_edgecolor(expression_values.color[i])
        slices[i].set_linewidth(0.1)

    # add black outline to cicle
    circle = Circle(center, radius, fill=False, edgecolor="k", linewidth=1)
    ax.add_patch(circle)


def _get_expression_table(adata, gene, color_map, max_expression):
    """helper function to extract expression table"""
    # get colormap
    cmap = get_cmap(color_map)

    # normalize item number values to colormap
    norm = Normalize(vmin=0, vmax=max_expression)

    if adata.n_obs > 1:

        # convert X to array
        if type(adata.X) == ndarray:
            pass
        else:
            adata.X = adata.X.toarray()

        if adata.n_obs > 1:

            # get expression values from adata
            if adata.n_vars == 1:
                expression_values = Series(adata.X.flatten()).value_counts().to_frame()
            else:
                iLoc = adata.var_names.tolist().index(gene)
                expression_values = (
                    Series(adata.X[:, iLoc].flatten()).value_counts().to_frame()
                )
            expression_values.columns = ["counts"]
            expression_values["value"] = expression_values.index.tolist()
            expression_values.sort_values(
                ascending=True, by="value", axis=0, inplace=True
            )
            expression_values.index = range(expression_values.shape[0])

            # colormap possible values = viridis, jet, spectral
            colors = []
            for i in expression_values.value:
                color = cmap(norm(i))
                colors.append(color)

            expression_values["color"] = colors

    elif adata.n_obs == 1 and adata.n_vars > 1:

        iLoc = adata.var_names.tolist().index(gene)
        expression_values = Series(adata.X[iLoc].flatten()).value_counts().to_frame()
        expression_values.columns = ["counts"]
        expression_values["value"] = expression_values.index.tolist()
        expression_values.sort_values(ascending=True, by="value", axis=0, inplace=True)
        expression_values.index = range(expression_values.shape[0])

        # colormap possible values = viridis, jet, spectral
        colors = []
        for i in expression_values.value:
            color = cmap(norm(i))
            colors.append(color)

        expression_values["color"] = colors

    elif adata.n_obs == 1 and adata.n_vars == 1:

        expression_values = DataFrame(
            {"counts": [1], "value": adata.X}, index=adata.obs.index
        )
        expression_values.sort_values(ascending=True, by="value", axis=0, inplace=True)
        expression_values.index = range(expression_values.shape[0])

        # colormap possible values = viridis, jet, spectral
        colors = []
        for i in expression_values.value:
            color = cmap(norm(i))
            colors.append(color)

        expression_values["color"] = colors

    else:
        # set to none since it doesnt contain any cells
        expression_values = None

    return expression_values


def _round_to_10(x):
    return 10 ** (ceil(log10(x)))


def _rescale(val, in_min, in_max, out_min, out_max):
    return out_min + (val - in_min) * ((out_max - out_min) / (in_max - in_min))


def dot_heatmap(
    adata,
    genes,
    group_by="louvain",
    plot_absolute=True,
    raw=True,
    color_map="Reds",
    figsize=None,
    order=None,
    switch_axis=False,
):
    """Generate a dot plot, filled with heatmap of individuals cells gene expression.

    This function generates a plot where the genes are plotted on the X-axis and the groups on
    the y-axis. For each coordinate a circle plot is generated where the size of the circle
    represents the number of cells in the group and the fill of the circle shows the gene
    expression of each individual cell as a sorted heatmap.

    If switch_axis = True then the genes are plotted on the y-axis and the groups on the x-axis.

    Note this function takes some time to calculate depending on the number of individual points
    that need to be plotted.

    parameters
    ----------
    adata: `AnnData`
        anndata object containing data that is to be visualized
    genes: `['str']`
        list of strings identifying the genes that are to be plotted
    group_by: `str` | default = 'louvain'
        string identifying the column in adata.obs that is to be used to generate the groups
    plot_absolute: `bool` | default
        boolian value used to indicate if the absolute geneexpression or the relativ gene expression
        should be plotted (currently not fully implemented)
    raw: `bool` | default = True
        boolian value indicating if the data saved in adata.raw should be used to generate the plot
    color_map: `str`
        string identifying the type of color_map that should be used to plot the results (needs to exist!)
    figsize: (value, value) or None | default = None
        optional parameter to define the figure size of the plot that is to be generated
    order: `['str']` or None | default = None
        optional parameter to pass an ordered list of groups to the function to specify the order in which
        they are to be plotted
    switch_axis: `bool` | default = False
        boolian value to determine if the axes should be switched. Per default the genes are plotted on the
        x-axis.

    returns
    -------
    Figure
        A matplotlib figure element containing the generated plot. To save the figure this plot will need
        to be passed to a parameter and saved in a second step through the fig.savefig() function call.

    Examples
    --------
    >>> pytest.skip('Test is only for here as example and should not be executed')
    >>> # import libraries and dataset
    >>> import besca as bc
    >>> adata = bc.datasets.Kotliarov2020_processed()
    >>> genes = ['CD3D', 'CD19']
    >>> fig = bc.pl.dot_heatmap(adata, genes=genes, group_by='leiden')

    .. plot::
        >>> pytest.skip('Test is only for here as example and should not be executed')
        >>> # import libraries and dataset
        >>> import besca as bc
        >>> adata = bc.datasets.Kotliarov2020_processed()
        >>> # define genes
        >>> genes = ['CD3D', 'CD19']
        >>> fig = bc.pl.dot_heatmap(adata, genes=genes, group_by='leiden')

    """
    # check to ensure genes does not contain any duplicates
    # if it does return warning message
    genes_orig = genes
    genes = [ele for ind, ele in enumerate(genes) if ele not in genes[:ind] or not ele]

    if len(genes) != len(genes_orig):
        print("WARNING: list of passed genes contained duplicates! These were removed.")

        # get duplicates
        seen = set()
        duplicates = []
        for ele in genes_orig:
            if ele in seen:
                duplicates.append(ele)
                break
            if ele:
                seen.add(ele)

        print("Duplicate genes: ", duplicates)

    if order is None:
        groups = adata.obs.get(group_by).value_counts().index.tolist()
    else:
        groups = order

    # reorder louvain groups to be in correct order
    if group_by == "louvain" and order is None:
        groups = [int(y) for y in groups]
        groups.sort()
        groups = [str(y) for y in groups]

    # get data to plot
    if raw:
        adata_plot = get_raw(adata)
        adata_plot.var_names_make_unique()
        gene_indexes = [adata_plot.var_names.tolist().index(x) for x in genes]
        adata_plot = adata_plot[:, gene_indexes]
    else:
        adata_plot = adata.copy()
        adata_plot.var_names_make_unique()
        gene_indexes = [adata_plot.var_names.tolist().index(x) for x in genes]
        adata_plot = adata_plot[:, gene_indexes]
        # get groups to be plotted on y-axis

    # get number of columns/rows
    if switch_axis:
        n_rows = len(genes) + 1  # additional rows: 1 row for the gene labels
        n_cols = (
            len(groups) + 6
        )  # additional columns: 2 columns as spacers, 1 column for group labels, 1 column for colorbar, 1 coloum for dotsize legend
    else:
        n_cols = (
            len(genes) + 6
        )  # additional columns: 2 columns as spacers, 1 column for group labels, 1 column for colorbar, 1 coloum for dotsize legend
        n_rows = len(groups) + 1  # additional rows: 1 row for the gene labels

    ############################
    # calculate bins for our data according to cellnumber in each group
    # get total cells
    total_num_cells = adata_plot.n_obs

    # get number of cells in each group
    n_cells = adata_plot.obs.get(group_by).value_counts().to_frame()
    n_cells.columns = ["n_cells"]

    # bin each group into a category
    n_cells["bin"] = [_round_to_10(x) for x in n_cells.n_cells]
    n_cells.bin = n_cells.bin.astype("category")

    bins = n_cells.bin.value_counts().index.tolist()

    if 1.0 not in bins:
        bins.append(1.0)
        bins.sort(reverse=True)

    n_bins = len(bins)
    range_per_bin = round(1.0 / n_bins, 6)

    value = 1.0
    values = [value]
    while value >= 0:
        value = round(value - range_per_bin, 5)
        if value >= 0:
            values.append(value)
    values = [round(x, 2) for x in values]

    # calculate properties of our bins
    bins = DataFrame(data={"bins": bins}, index=bins)
    bins.sort_values(ascending=False, axis=0, by="bins", inplace=True)
    bins["upper_limit_in"] = bins.bins
    bins["lower_limit_in"] = bins.bins.tolist()[1:] + [-1]
    bins["lower_limit_in"] = (bins.lower_limit_in + 1).tolist()

    bins["upper_limit_out"] = values[:-1]
    bins["lower_limit_out"] = bins.upper_limit_out.tolist()[1:] + [-0.01]
    bins["lower_limit_out"] = (bins.lower_limit_out + 0.01).tolist()

    n_cells["radius"] = [
        _rescale(
            val=x,
            in_max=bins["upper_limit_in"].get(bin_id),
            in_min=bins["lower_limit_in"].get(bin_id),
            out_max=bins["upper_limit_out"].get(bin_id),
            out_min=bins["lower_limit_out"].get(bin_id),
        )
        for x, bin_id in zip(n_cells.n_cells.tolist(), n_cells.bin.tolist())
    ]

    # check that we dont have more bins than rows
    if n_bins + 1 > n_rows:
        if switch_axis:
            while len(genes) - n_bins < 0:
                genes.append("NA")
        else:
            while len(groups) - n_bins < 0:
                groups.append("NA")
        n_rows = n_bins + 1

    # setup plot dimensions
    colorbar_width = 0.2
    colorbar_width_spacer = 0.5
    size_legend_spacer = 1
    size_legend_width = 0.5
    size_legend_label = 0.7
    label_width = 3
    label_height = 3

    if figsize is None:
        if switch_axis:
            heatmap_height = len(genes) * 0.5
            heatmap_width = len(groups) * 0.5
        else:
            heatmap_height = len(groups) * 0.5
            heatmap_width = len(genes) * 0.5
        heatmap_height = max(
            [1.5, heatmap_height]
        )  # ensure that we have a certain minimal height even for small number of groups
        height = heatmap_height + label_height
        width = (
            label_width
            + heatmap_width
            + colorbar_width_spacer
            + colorbar_width
            + size_legend_spacer
            + size_legend_width
            + size_legend_label
        )
    else:
        width, height = figsize
        heatmap_width = width - (
            colorbar_width_spacer
            + colorbar_width
            + size_legend_spacer
            + size_legend_width
            + size_legend_label
        )
        heatmap_heigth = height - label_height

    if switch_axis:
        width_group = heatmap_width / len(groups)
        height_gene = heatmap_height / len(genes)
    else:
        width_gene = heatmap_width / len(genes)
        height_group = heatmap_height / len(groups)

    # determine ratios
    if switch_axis:
        width_ratios = (
            [label_width]
            + [width_group] * len(groups)
            + [
                colorbar_width_spacer,
                colorbar_width,
                size_legend_spacer,
                size_legend_width,
                size_legend_label,
            ]
        )
        height_ratios = [label_height] + [height_gene] * len(genes)
    else:
        width_ratios = (
            [label_width]
            + [width_gene] * len(genes)
            + [
                colorbar_width_spacer,
                colorbar_width,
                size_legend_spacer,
                size_legend_width,
                size_legend_label,
            ]
        )
        height_ratios = [label_height] + [height_group] * len(groups)

    ##### generate the base layout of our plot
    # nrows = 1 for labeling, + number of groups that are to be plotted
    # ncols = 1 for labeling, + number of genes that are to be plotted, + 1 for colorbar, +1 as spacer, +1 as size legend
    fig = figure(figsize=(width, height))
    axs = gridspec.GridSpec(
        nrows=n_rows,
        ncols=n_cols,
        wspace=0.02,
        hspace=0.04,
        width_ratios=width_ratios,
        height_ratios=height_ratios,
    )

    # adjust height of the colorlegend so it does not go the entire length of the figure for a large number of rows
    if n_rows < 8:
        color_legend = fig.add_subplot(axs[1:, n_cols - 4])
    else:
        colorbar_height = min(5.0, height - label_height)
        wspace = 10.5 / width
        axs2 = gridspec.GridSpecFromSubplotSpec(
            2,
            1,
            subplot_spec=axs[1:, n_cols - 4],
            wspace=wspace,
            height_ratios=[
                colorbar_height / height,
                (height - colorbar_height) / height,
            ],
        )

        # add ColorBarlegend
        color_legend = fig.add_subplot(axs2[0])

    if plot_absolute:
        max_expression = adata_plot.X.max()
    else:
        max_expression = 1

    cmap = get_cmap(color_map)
    norm = Normalize(vmin=0, vmax=max_expression)
    ColorbarBase(color_legend, cmap=cmap, norm=norm)

    # add legend name
    label_num_cells = fig.add_subplot(axs[0, -4])
    label_num_cells.axis("off")
    text = label_num_cells.text(
        0.2, 0.1, "log2cp10k", ha="center", va="bottom", rotation=0
    )

    # add dotsize legend relating to cellnumber in group
    label_num_cells = fig.add_subplot(axs[0, -2])
    label_num_cells.axis("off")
    text = label_num_cells.text(
        0.2, 0.1, "number of cells", ha="left", va="bottom", rotation=0
    )

    for n in range(n_bins):
        this_bin = bins.bins.tolist()[n]
        radius = bins.upper_limit_out.get(this_bin)

        bin_plot = fig.add_subplot(axs[n + 1, -2])
        bin_plot.axis("equal")
        sizes = [1]
        bin_plot.pie(
            sizes,
            colors=["white"],
            shadow=False,
            center=(0.5, 0.5),
            radius=radius,
            wedgeprops={"edgecolor": "k", "linewidth": 1},
        )

        bin_label = fig.add_subplot(axs[n + 1, -1])
        bin_label.axis("off")
        text = bin_label.text(
            0.2,
            0.5,
            str(int(bins.upper_limit_in.get(this_bin))),
            ha="left",
            va="center",
            rotation=0,
        )

    # add expression data

    genes = [value for value in genes if value != "NA"]
    groups = [value for value in groups if value != "NA"]

    if switch_axis:
        sub_num_rows = len(genes)
        sub_num_cols = len(groups)
    else:
        sub_num_rows = len(groups)
        sub_num_cols = len(genes)

    if switch_axis:
        for i in range(sub_num_cols):
            adata_subset = adata_plot[
                adata_plot.obs.get(group_by) == groups[i], :
            ].copy()
            adata_subset.var_names_make_unique()
            label = fig.add_subplot(axs[0, i + 1])
            label.axis("off")
            text = label.text(
                0.5, 0.2, groups[i], ha="center", va="bottom", rotation=90
            )
            for j in range(sub_num_rows):
                data = _get_expression_table(
                    adata_subset,
                    genes[j],
                    color_map=color_map,
                    max_expression=max_expression,
                )
                if data is not None:
                    radius = n_cells.radius.get(groups[i])
                    ax1 = fig.add_subplot(axs[j + 1, i + 1])
                    ax1.axis("equal")
                    _generate_circle(data, (0.5, 0.5), radius, ax=ax1)
    else:
        for i in range(sub_num_rows):
            adata_subset = adata_plot[
                adata_plot.obs.get(group_by) == groups[i], :
            ].copy()
            adata_subset.var_names_make_unique()
            label = fig.add_subplot(axs[i + 1, 0])
            label.axis("off")
            text = label.text(0.8, 0.5, groups[i], ha="right", va="center")
            for j in range(sub_num_cols):
                data = _get_expression_table(
                    adata_subset,
                    genes[j],
                    color_map=color_map,
                    max_expression=max_expression,
                )
                if data is not None:
                    radius = n_cells.radius.get(groups[i])
                    ax1 = fig.add_subplot(axs[i + 1, j + 1])
                    ax1.axis("equal")
                    _generate_circle(data, (0.5, 0.5), radius, ax=ax1)

    # add gene labels
    if switch_axis:
        for j in range(sub_num_rows):
            label = fig.add_subplot(axs[j + 1, 0])
            label.axis("off")
            text = label.text(0.8, 0.5, genes[j], ha="right", va="center")
    else:
        for j in range(sub_num_cols):
            label = fig.add_subplot(axs[0, j + 1])
            label.axis("off")
            text = label.text(0.5, 0.2, genes[j], ha="center", va="bottom", rotation=90)

    fig.tight_layout()
    return fig


def dot_heatmap_split(
    adata,
    genes,
    split_by,
    group_by="leiden",
    plot_absolute=True,
    raw=True,
    color_maps=["Reds", "Blues"],
    figsize=None,
    order=None,
    switch_axis=False,
    rotation=None,
):
    """Generate a dot plot, filled with heatmap of individuals cells gene expression to compare two conditions.

    This function generates a plot where the genes are plotted on the X-axis and the groups on
    the y-axis. For each coordinate two side-by-side circle plots are generated where the size of the circle
    represents the number of cells in the group and the fill of the circle shows the gene
    expression of each individual cell as a sorted heatmap. The two circles represent the cells in that group
    from each condition.

    If switch_axis = True then the genes are plotted on the y-axis and the groups on the x-axis.

    Note this function takes some time to calculate depending on the number of individual points
    that need to be plotted.

    Note: currently this function is only implemented to support pairwise comparisions (i.e. with a maximum of
    two groups)

    parameters
    ----------
    adata: `AnnData`
        anndata object containing data that is to be visualized
    genes: `['str']`
        list of strings identifying the genes that are to be plotted
    split_by: `str`
        string identifying the column in adata.obs that is to be used to identify the two conditions.
        May only contain two different labels.
    group_by: `str` | default = 'louvain'
        string identifying the column in adata.obs that is to be used to generate the groups
    plot_absolute: `bool` | default
        boolian value used to indicate if the absolute geneexpression or the relativ gene expression
        should be plotted (currently not fully implemented)
    raw: `bool` | default = True
        boolian value indicating if the data saved in adata.raw should be used to generate the plot
    color_map: `str`
        string identifying the type of color_map that should be used to plot the results (needs to exist!)
    figsize: (value, value) or None | default = None
        optional parameter to define the figure size of the plot that is to be generated
    order: `['str']` or None | default = None
        optional parameter to pass an ordered list of groups to the function to specify the order in which
        they are to be plotted
    switch_axis: `bool` | default = False
        boolian value to determine if the axes should be switched. Per default the genes are plotted on the
        x-axis.
    rotation: `int`
        optional parameter to define by which degree the labels on the x-axis should be rotatet. Per default
        this is 0 if the genes are plotted on the x-axis and 90 if the group labels are plotted on the x-axis.

    returns
    -------
    Figure
        A matplotlib figure element containing the generated plot. To save the figure this plot will need
        to be passed to a parameter and saved in a second step through the fig.savefig() function call.

    Examples
    --------
    >>> pytest.skip('Test is only for here as example and should not be executed')
    >>> # import libraries and dataset
    >>> import besca as bc
    >>> adata = bc.datasets.Haber2017_processed()
    >>> genes = ['Defa22', 'Defa24', 'Gm15284', 'Reg4']
    >>> fig = bc.pl.dot_heatmap_split(adata, genes=genes, group_by='leiden', split_by = 'donor')

    .. plot::
        >>> pytest.skip('Test is only for here as example and should not be executed')
        >>> # import libraries and dataset
        >>> import besca as bc
        >>> adata = bc.datasets.Haber2017_processed()
        >>> # define genes
        >>> genes = ['Defa22', 'Defa24', 'Gm15284', 'Reg4']
        >>> fig = bc.pl.dot_heatmap_split(adata, genes=genes, group_by='leiden', split_by = 'donor')

    """

    # check to ensure genes does not contain any duplicates
    # if it does return warning message
    genes_orig = genes
    genes = [ele for ind, ele in enumerate(genes) if ele not in genes[:ind] or not ele]

    if len(genes) != len(genes_orig):
        print("WARNING: list of passed genes contained duplicates! These were removed.")

        # get duplicates
        seen = set()
        duplicates = []
        for ele in genes_orig:
            if ele in seen:
                duplicates.append(ele)
                break
            if ele:
                seen.add(ele)

        print("Duplicate genes: ", duplicates)

    if order is None:
        groups = adata.obs.get(group_by).value_counts().index.tolist()
    else:
        groups = order

    # reorder louvain groups to be in correct order
    if group_by == "louvain" and order is None:
        groups = [int(y) for y in groups]
        groups.sort()
        groups = [str(y) for y in groups]

    # get data to plot
    if raw:
        adata_plot = get_raw(adata)
        adata_plot.var_names_make_unique()
        gene_indexes = [adata_plot.var_names.tolist().index(x) for x in genes]
        adata_plot = adata_plot[:, gene_indexes].copy()
    else:
        adata_plot = adata.copy()
        adata_plot.var_names_make_unique()
        gene_indexes = [adata_plot.var_names.tolist().index(x) for x in genes]
        adata_plot = adata_plot[:, gene_indexes].copy()
        # get groups to be plotted on y-axis

    # get data for the groups we want to plot against each other
    if adata.obs.get(split_by).dtype.name == "category":
        conditions = adata.obs.get(split_by).cat.categories.tolist()
    else:
        conditions = adata_plot.obs.get(split_by).value_counts().index.tolist()

    n_conditions = len(conditions)

    if n_conditions != 2:
        sys.exit(
            "Can only compare two conditions! Please ensure that adata.obs.get(split_variable) only contains two conditions."
        )

    # generate data subsets
    condition1 = adata_plot[adata_plot.obs.get(split_by) == conditions[0], :]
    condition1.var_names_make_unique()
    condition2 = adata_plot[adata_plot.obs.get(split_by) == conditions[1], :]
    condition2.var_names_make_unique()

    # get number of columns/rows
    if switch_axis:
        if rotation is None:
            rotation = 90
        n_rows = len(genes) + 1  # additional rows: 1 row for the gene labels
        n_cols = (
            2 * len(groups) + len(groups) - 1 + 8
        )  # additional columns:(groups - 1) X dis_columns between groups, 2 columns as spacers, 1 column for group labels, 1 column for colorbar, 1 coloum for dotsize legend
    else:
        if rotation is None:
            rotation = 0
        n_cols = (
            2 * len(genes) + len(genes) - 1 + 8
        )  # additional columns: (genes - 1) X dis_columns between groups, 2 columns as spacers, 1 column for group labels, 1 column for colorbar, 1 coloum for dotsize legend
        n_rows = len(groups) + 1  # additional rows: 1 row for the gene labels

    ############################
    # calculate bins for our data according to cellnumber in each group

    # get number of cells in each group
    n_cells_cond1 = condition1.obs.get(group_by).value_counts().to_frame()
    n_cells_cond1.columns = ["n_cells"]

    n_cells_cond2 = condition2.obs.get(group_by).value_counts().to_frame()
    n_cells_cond2.columns = ["n_cells"]

    # bin each group into a category
    n_cells_cond1["bin"] = [_round_to_10(x) for x in n_cells_cond1.n_cells]
    n_cells_cond1.bin = n_cells_cond1.bin.astype("category")

    n_cells_cond2["bin"] = [_round_to_10(x) for x in n_cells_cond2.n_cells]
    n_cells_cond2.bin = n_cells_cond2.bin.astype("category")

    # get all bins from both datasets
    bins = list(
        set(
            n_cells_cond1.bin.value_counts().index.tolist()
            + n_cells_cond2.bin.value_counts().index.tolist()
        )
    )

    # ensure that smallest circle always represents 1 cell (ensures that scales for all circles are large enough)
    if 1.0 not in bins:
        bins.append(1.0)
        bins.sort(reverse=True)

    n_bins = len(bins)
    range_per_bin = round(1.0 / n_bins, 6)

    value = 1.0
    values = [value]
    while value >= 0:
        value = round(value - range_per_bin, 5)
        if value >= 0:
            values.append(value)
    values = [round(x, 2) for x in values]

    # calculate properties of our bins
    bins = DataFrame(data={"bins": bins}, index=bins)
    bins.sort_values(ascending=False, axis=0, by="bins", inplace=True)
    bins["upper_limit_in"] = bins.bins
    bins["lower_limit_in"] = bins.bins.tolist()[1:] + [-1]
    bins["lower_limit_in"] = (bins.lower_limit_in + 1).tolist()

    bins["upper_limit_out"] = values[:-1]
    bins["lower_limit_out"] = bins.upper_limit_out.tolist()[1:] + [-0.01]
    bins["lower_limit_out"] = (bins.lower_limit_out + 0.01).tolist()

    # calculate radius for each group
    n_cells_cond1["radius"] = [
        _rescale(
            val=x,
            in_max=bins["upper_limit_in"].get(bin_id),
            in_min=bins["lower_limit_in"].get(bin_id),
            out_max=bins["upper_limit_out"].get(bin_id),
            out_min=bins["lower_limit_out"].get(bin_id),
        )
        for x, bin_id in zip(n_cells_cond1.n_cells.tolist(), n_cells_cond1.bin.tolist())
    ]

    n_cells_cond2["radius"] = [
        _rescale(
            val=x,
            in_max=bins["upper_limit_in"].get(bin_id),
            in_min=bins["lower_limit_in"].get(bin_id),
            out_max=bins["upper_limit_out"].get(bin_id),
            out_min=bins["lower_limit_out"].get(bin_id),
        )
        for x, bin_id in zip(n_cells_cond2.n_cells.tolist(), n_cells_cond2.bin.tolist())
    ]

    # check that we dont have more bins than rows
    # if we do momentarily append 'NA's to genes to get correct number of rows
    # 'NA's removed later

    if n_bins + 1 > n_rows:
        if switch_axis:
            while len(genes) - n_bins < 0:
                genes.append("NA")
        else:
            while len(groups) - n_bins < 0:
                groups.append("NA")
        n_rows = n_bins + 1

    # setup plot dimensions
    colorbar_width = 0.2
    colorbar_width_spacer = 0.5
    between_colorbar_spacer = 0.8
    size_legend_spacer = 1
    size_legend_width = 0.5
    size_legend_label = 0.7
    label_width = 3
    label_height = 3
    dis_pairs = 0.1

    if figsize is None:
        if switch_axis:
            heatmap_height = len(genes) * 0.5
            heatmap_width = 2 * len(groups) * 0.5 + dis_pairs * (len(groups) - 1)
        else:
            heatmap_height = len(groups) * 0.5
            heatmap_width = 2 * len(genes) * 0.5 + dis_pairs * (len(genes) - 1)
        heatmap_height = max(
            [1.5, heatmap_height]
        )  # ensure that we have a certain minimal height even for small number of groups
        height = heatmap_height + label_height
        width = (
            label_width
            + heatmap_width
            + colorbar_width_spacer
            + colorbar_width
            + between_colorbar_spacer
            + colorbar_width
            + size_legend_spacer
            + size_legend_width
            + size_legend_label
        )
    else:
        width, height = figsize
        heatmap_width = width - (
            label_width
            + colorbar_width_spacer
            + colorbar_width
            + between_colorbar_spacer
            + colorbar_width
            + size_legend_spacer
            + size_legend_width
            + size_legend_label
        )
        heatmap_heigth = height - label_height

    if switch_axis:
        width_group = (heatmap_width - dis_pairs * (len(groups) - 1)) / (
            2 * len(groups)
        )
        height_gene = heatmap_height / len(genes)
    else:
        width_gene = (heatmap_width - dis_pairs * (len(genes) - 1)) / (2 * len(genes))
        height_group = heatmap_height / len(groups)

    # determine ratios
    if switch_axis:
        width_ratios = (
            [label_width]
            + [width_group] * 2
            + [dis_pairs, width_group, width_group] * (len(groups) - 1)
            + [
                colorbar_width_spacer,
                colorbar_width,
                between_colorbar_spacer,
                colorbar_width,
                size_legend_spacer,
                size_legend_width,
                size_legend_label,
            ]
        )
        height_ratios = [label_height] + [height_gene] * len(genes)
    else:
        width_ratios = (
            [label_width]
            + [width_gene] * 2
            + [dis_pairs, width_gene, width_gene] * (len(genes) - 1)
            + [
                colorbar_width_spacer,
                colorbar_width,
                between_colorbar_spacer,
                colorbar_width,
                size_legend_spacer,
                size_legend_width,
                size_legend_label,
            ]
        )
        height_ratios = [label_height] + [height_group] * len(groups)

    ##### generate the base layout of our plot
    # nrows = 1 for labeling, + number of groups that are to be plotted
    # ncols = 1 for labeling, + number of genes that are to be plotted, + 1 for colorbar, +1 as spacer, +1 as size legend
    fig = figure(figsize=(width, height))
    axs = gridspec.GridSpec(
        nrows=n_rows,
        ncols=n_cols,
        wspace=0.02,
        hspace=0.04,
        width_ratios=width_ratios,
        height_ratios=height_ratios,
    )
    # adjust height of the colorlegend so it does not go the entire length of the figure for a large number of rows
    if n_rows < 8:
        color_legend_cond1 = fig.add_subplot(axs[2:, -6])
        color_legend_cond2 = fig.add_subplot(axs[2:, -4])
    else:
        colorbar_height = min(5.0, height - label_height)
        wspace = 10.5 / width
        axs2_cond1 = gridspec.GridSpecFromSubplotSpec(
            2,
            1,
            subplot_spec=axs[2:, -6],
            wspace=wspace,
            height_ratios=[
                colorbar_height / height,
                (height - colorbar_height) / height,
            ],
        )
        axs2_cond2 = gridspec.GridSpecFromSubplotSpec(
            2,
            1,
            subplot_spec=axs[2:, -4],
            wspace=wspace,
            height_ratios=[
                colorbar_height / height,
                (height - colorbar_height) / height,
            ],
        )

        # add ColorBarlegend
        color_legend_cond1 = fig.add_subplot(axs2_cond1[0])
        color_legend_cond2 = fig.add_subplot(axs2_cond2[0])

    if plot_absolute:
        max_expression = adata_plot.X.max()
    else:
        max_expression = 1

    cmap_cond1 = get_cmap(color_maps[0])
    cmap_cond2 = get_cmap(color_maps[1])
    norm = Normalize(vmin=0, vmax=max_expression)
    ColorbarBase(color_legend_cond1, cmap=cmap_cond1, norm=norm)
    ColorbarBase(color_legend_cond2, cmap=cmap_cond2, norm=norm)

    # add legend name (log2cp10k)
    label = fig.add_subplot(axs[0, -5])
    label.axis("off")
    text = label.text(
        0.5, 0.05, "log2cp10k", ha="center", va="bottom", rotation=0, fontsize="x-large"
    )

    # add legend for each condition
    label_cond1 = fig.add_subplot(axs[1, -6])
    label_cond1.axis("off")
    text = label_cond1.text(
        0.5, 0.5, conditions[0], ha="center", va="center", rotation=0
    )

    label_cond2 = fig.add_subplot(axs[1, -4])
    label_cond2.axis("off")
    text = label_cond2.text(
        0.5, 0.5, conditions[1], ha="center", va="center", rotation=0
    )

    # add dotsize legend relating to cellnumber in group
    label_num_cells = fig.add_subplot(axs[0, -2])
    label_num_cells.axis("off")
    text = label_num_cells.text(
        0.0,
        0.05,
        "number of cells",
        ha="left",
        va="bottom",
        rotation=0,
        fontsize="x-large",
    )

    for n in range(n_bins):
        this_bin = bins.bins.tolist()[n]
        radius = bins.upper_limit_out.get(this_bin)

        bin_plot = fig.add_subplot(axs[n + 1, -2])
        bin_plot.axis("equal")
        sizes = [1]
        bin_plot.pie(
            sizes,
            colors=["white"],
            shadow=False,
            center=(0.5, 0.5),
            radius=radius,
            wedgeprops={"edgecolor": "k", "linewidth": 1},
        )

        bin_label = fig.add_subplot(axs[n + 1, -1])
        bin_label.axis("off")
        text = bin_label.text(
            0.2,
            0.5,
            str(int(bins.upper_limit_in.get(this_bin))),
            ha="left",
            va="center",
            rotation=0,
        )

    # add expression data
    genes = [value for value in genes if value != "NA"]
    groups = [value for value in groups if value != "NA"]

    if switch_axis:
        sub_num_rows = len(genes)
        sub_num_cols = 2 * len(groups) + len(groups) - 1
    else:
        sub_num_rows = len(groups)
        sub_num_cols = 2 * len(genes) + len(genes) - 1

    col_range_cond1 = list(arange(1, sub_num_cols + 1, 3))
    col_range_cond2 = list(arange(2, sub_num_cols + 1, 3))
    colranges = [col_range_cond1, col_range_cond2]

    input_data = {conditions[0]: condition1, conditions[1]: condition2}
    n_cells = {conditions[0]: n_cells_cond1, conditions[1]: n_cells_cond2}

    if switch_axis:
        for c in range(len(conditions)):
            # get datasubset containing only that condition
            adata_subset_plot = input_data[conditions[c]]
            # iterate through the groups (here on the x-axis)
            for i, g in zip(colranges[c], range(len(groups))):
                adata_subset = adata_subset_plot[
                    adata_subset_plot.obs.get(group_by) == groups[g], :
                ].copy()
                # iterate through all of the genes
                for j in range(sub_num_rows):
                    data = _get_expression_table(
                        adata_subset,
                        genes[j],
                        color_map=color_maps[c],
                        max_expression=max_expression,
                    )
                    if data is not None:
                        radius = n_cells[conditions[c]].radius.get(groups[g])
                        ax1 = fig.add_subplot(axs[j + 1, i])
                        ax1.axis("equal")
                        _generate_circle(data, (0.5, 0.5), radius, ax=ax1)
    else:
        for c in range(len(conditions)):
            # get datasubset containing only that condition
            adata_subset_plot = input_data[conditions[c]]
            # iterate through the groups (here on the y-axis)
            for i in range(sub_num_rows):
                adata_subset = adata_subset_plot[
                    adata_subset_plot.obs.get(group_by) == groups[i], :
                ].copy()
                # iterate through all of the genes
                for j, g in zip(colranges[c], range(len(genes))):
                    data = _get_expression_table(
                        adata_subset,
                        genes[g],
                        color_map=color_maps[c],
                        max_expression=max_expression,
                    )
                    if data is not None:
                        radius = n_cells[conditions[c]].radius.get(groups[i])
                        ax1 = fig.add_subplot(axs[i + 1, j])
                        ax1.axis("equal")
                        _generate_circle(data, (0.5, 0.5), radius, ax=ax1)

    # add group labels
    if switch_axis:
        for j, g in zip(colranges[0], range(len(groups))):
            label = fig.add_subplot(axs[0, j : (j + 2)])
            label.axis("off")
            text = label.text(
                0.5, 0.05, groups[g], ha="center", va="bottom", rotation=rotation
            )
            label.annotate(
                "",
                xy=(0, 0),
                xycoords="data",
                xytext=(1, 0),
                textcoords="data",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=0."),
            )
    else:
        for j in range(sub_num_rows):
            label = fig.add_subplot(axs[j + 1, 0])
            label.axis("off")
            text = label.text(0.95, 0.5, groups[j], ha="right", va="center")

    # add gene labels
    if switch_axis:
        for j in range(sub_num_rows):
            label = fig.add_subplot(axs[j + 1, 0])
            label.axis("off")
            text = label.text(0.95, 0.5, genes[j], ha="right", va="center")
    else:
        for j, g in zip(colranges[0], range(len(genes))):
            label = fig.add_subplot(axs[0, j : (j + 2)])
            label.axis("off")
            text = label.text(
                0.5, 0.05, genes[g], ha="center", va="bottom", rotation=rotation
            )
            label.annotate(
                "",
                xy=(0, 0),
                xycoords="data",
                xytext=(1, 0),
                textcoords="data",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=0."),
            )
    fig.tight_layout()
    return fig


def dot_heatmap_split_greyscale(
    adata,
    genes,
    split_by,
    group_by="louvain",
    abbreviations=None,
    plot_absolute=True,
    raw=True,
    figsize=None,
    order=None,
    switch_axis=False,
    rotation=None,
):
    """Generate a dot plot, filled with heatmap of individuals cells gene expression to compare two conditions (greyscale).

    This function generates a plot where the genes are plotted on the X-axis and the groups on
    the y-axis. For each coordinate two side-by-side circle plots are generated where the size of the circle
    represents the number of cells in the group and the fill of the circle shows the gene
    expression of each individual cell as a sorted heatmap. The two circles represent the cells in that group
    from each condition.

    If switch_axis = True then the genes are plotted on the y-axis and the groups on the x-axis.

    Note this function takes some time to calculate depending on the number of individual points
    that need to be plotted.

    Note: currently this function is only implemented to support pairwise comparisions (i.e. with a maximum of
    two groups)

    parameters
    ----------
    adata: `AnnData`
        anndata object containing data that is to be visualized
    genes: `['str']`
        list of strings identifying the genes that are to be plotted
    split_by: `str`
        string identifying the column in adata.obs that is to be used to identify the two conditions.
        May only contain two different labels.
    group_by: `str` | default = 'louvain'
        string identifying the column in adata.obs that is to be used to generate the groups
    abbreviations: `['str']` | default = None
        optional parameter which accepts a list of strings identifying the abbreviations for the two groups that
        should be plotted atop each dot row, if none is supplied it plots the entire group name. The list needs
        to be supplied in the correct order assuming alphabetically sorted conditions.
    plot_absolute: `bool` | default
        boolian value used to indicate if the absolute geneexpression or the relativ gene expression
        should be plotted (currently not fully implemented)
    raw: `bool` | default = True
        boolian value indicating if the data saved in adata.raw should be used to generate the plot
    figsize: (value, value) or None | default = None
        optional parameter to define the figure size of the plot that is to be generated
    order: `['str']` or None | default = None
        optional parameter to pass an ordered list of groups to the function to specify the order in which
        they are to be plotted
    switch_axis: `bool` | default = False
        boolian value to determine if the axes should be switched. Per default the genes are plotted on the
        x-axis.
    rotation: `int`
        optional parameter to define by which degree the labels on the x-axis should be rotatet. Per default
        this is 0 if the genes are plotted on the x-axis and 90 if the group labels are plotted on the x-axis.

    returns
    -------
    Figure
        A matplotlib figure element containing the generated plot. To save the figure this plot will need
        to be passed to a parameter and saved in a second step through the fig.savefig() function call.


    """

    # check to ensure genes does not contain any duplicates
    # if it does return warning message
    genes_orig = genes
    genes = [ele for ind, ele in enumerate(genes) if ele not in genes[:ind] or not ele]

    if len(genes) != len(genes_orig):
        print("WARNING: list of passed genes contained duplicates! These were removed.")

        # get duplicates
        seen = set()
        duplicates = []
        for ele in genes_orig:
            if ele in seen:
                duplicates.append(ele)
                break
            if ele:
                seen.add(ele)

        print("Duplicate genes: ", duplicates)

    if order is None:
        groups = adata.obs.get(group_by).value_counts().index.tolist()
    else:
        groups = order

    # reorder louvain groups to be in correct order
    if group_by == "louvain" and order is None:
        groups = [int(y) for y in groups]
        groups.sort()
        groups = [str(y) for y in groups]

    # get data to plot
    if raw:
        adata_plot = get_raw(adata)
        adata_plot.var_names_make_unique()
        gene_indexes = [adata_plot.var_names.tolist().index(x) for x in genes]
        adata_plot = adata_plot[:, gene_indexes].copy()
    else:
        adata_plot = adata.copy()
        adata_plot.var_names_make_unique()
        gene_indexes = [adata_plot.var_names.tolist().index(x) for x in genes]
        adata_plot = adata_plot[:, gene_indexes].copy()
        # get groups to be plotted on y-axis

    # get data for the groups we want to plot against each other
    conditions = adata_plot.obs.get(split_by).value_counts().index.tolist()
    # sort so they are alphabetic
    conditions.sort()

    n_conditions = len(conditions)

    if n_conditions != 2:
        sys.exit(
            "Can only compare two conditions! Please ensure that adata.obs.get(split_variable) only contains two conditions."
        )

    if abbreviations is None:
        abbrev_labels = conditions
    else:
        abbrev_labels = abbreviations

    # generate data subsets
    condition1 = adata_plot[adata_plot.obs.get(split_by) == conditions[0], :]
    condition1.var_names_make_unique()
    condition2 = adata_plot[adata_plot.obs.get(split_by) == conditions[1], :]
    condition2.var_names_make_unique()

    # get number of columns/rows
    if switch_axis:
        if rotation is None:
            rotation = 90
        n_rows = (
            len(genes) + 2
        )  # additional rows: 1 row for the gene labels and 1 for condition labels
        n_cols = (
            2 * len(groups) + len(groups) - 1 + 6
        )  # additional columns:(groups - 1) X dis_columns between groups, 2 columns as spacers, 1 column for group labels, 1 column for colorbar, 2 coloum for dotsize legend
    else:
        if rotation is None:
            rotation = 0
        n_cols = (
            2 * len(genes) + len(genes) - 1 + 6
        )  # additional columns: (genes - 1) X dis_columns between groups, 2 columns as spacers, 1 column for group labels, 1 column for colorbar, 1 coloum for dotsize legend
        n_rows = len(groups) + 2  # additional rows: 1 row for the gene labels

    ############################
    # calculate bins for our data according to cellnumber in each group

    # get number of cells in each group
    n_cells_cond1 = condition1.obs.get(group_by).value_counts().to_frame()
    n_cells_cond1.columns = ["n_cells"]

    n_cells_cond2 = condition2.obs.get(group_by).value_counts().to_frame()
    n_cells_cond2.columns = ["n_cells"]

    # bin each group into a category
    n_cells_cond1["bin"] = [_round_to_10(x) for x in n_cells_cond1.n_cells]
    n_cells_cond1.bin = n_cells_cond1.bin.astype("category")

    n_cells_cond2["bin"] = [_round_to_10(x) for x in n_cells_cond2.n_cells]
    n_cells_cond2.bin = n_cells_cond2.bin.astype("category")

    # get all bins from both datasets
    bins = list(
        set(
            n_cells_cond1.bin.value_counts().index.tolist()
            + n_cells_cond2.bin.value_counts().index.tolist()
        )
    )

    # ensure that smallest circle always represents 1 cell (ensures that scales for all circles are large enough)
    if 1.0 not in bins:
        bins.append(1.0)
        bins.sort(reverse=True)

    n_bins = len(bins)
    range_per_bin = round(1.0 / n_bins, 6)

    value = 1.0
    values = [value]
    while value >= 0:
        value = round(value - range_per_bin, 5)
        if value >= 0:
            values.append(value)
    values = [round(x, 2) for x in values]

    # calculate properties of our bins
    bins = DataFrame(data={"bins": bins}, index=bins)
    bins.sort_values(ascending=False, axis=0, by="bins", inplace=True)
    bins["upper_limit_in"] = bins.bins
    bins["lower_limit_in"] = bins.bins.tolist()[1:] + [-1]
    bins["lower_limit_in"] = (bins.lower_limit_in + 1).tolist()

    bins["upper_limit_out"] = values[:-1]
    bins["lower_limit_out"] = bins.upper_limit_out.tolist()[1:] + [-0.01]
    bins["lower_limit_out"] = (bins.lower_limit_out + 0.01).tolist()

    # calculate radius for each group
    n_cells_cond1["radius"] = [
        _rescale(
            val=x,
            in_max=bins["upper_limit_in"].get(bin_id),
            in_min=bins["lower_limit_in"].get(bin_id),
            out_max=bins["upper_limit_out"].get(bin_id),
            out_min=bins["lower_limit_out"].get(bin_id),
        )
        for x, bin_id in zip(n_cells_cond1.n_cells.tolist(), n_cells_cond1.bin.tolist())
    ]

    n_cells_cond2["radius"] = [
        _rescale(
            val=x,
            in_max=bins["upper_limit_in"].get(bin_id),
            in_min=bins["lower_limit_in"].get(bin_id),
            out_max=bins["upper_limit_out"].get(bin_id),
            out_min=bins["lower_limit_out"].get(bin_id),
        )
        for x, bin_id in zip(n_cells_cond2.n_cells.tolist(), n_cells_cond2.bin.tolist())
    ]

    # check that we dont have more bins than rows
    # if we do momentarily append 'NA's to genes to get correct number of rows
    # 'NA's removed later

    if n_bins + 1 > n_rows:
        if switch_axis:
            while len(genes) - n_bins < 0:
                genes.append("NA")
        else:
            while len(groups) - n_bins < 0:
                groups.append("NA")
        n_rows = n_bins + 1

    # setup plot dimensions
    colorbar_width = 0.2
    colorbar_width_spacer = 0.6
    size_legend_spacer = 0.7
    size_legend_width = 0.5
    size_legend_label = 0.7
    label_width = 3
    label_height = 3
    condition_label = 0.4
    dis_pairs = 0.1

    if figsize is None:
        if switch_axis:
            heatmap_height = len(genes) * 0.5
            heatmap_width = 2 * len(groups) * 0.5 + dis_pairs * (len(groups) - 1)
        else:
            heatmap_height = len(groups) * 0.5
            heatmap_width = 2 * len(genes) * 0.5 + dis_pairs * (len(genes) - 1)
        heatmap_height = max(
            [1.5, heatmap_height]
        )  # ensure that we have a certain minimal height even for small number of groups
        height = heatmap_height + label_height + condition_label
        width = (
            label_width
            + heatmap_width
            + colorbar_width_spacer
            + colorbar_width
            + size_legend_spacer
            + size_legend_width
            + size_legend_label
        )
    else:
        width, height = figsize
        heatmap_width = width - (
            label_width
            + colorbar_width_spacer
            + colorbar_width
            + size_legend_spacer
            + size_legend_width
            + size_legend_label
        )
        heatmap_heigth = height - label_height - condition_label

    if switch_axis:
        width_group = (heatmap_width - dis_pairs * (len(groups) - 1)) / (
            2 * len(groups)
        )
        height_gene = heatmap_height / len(genes)
    else:
        width_gene = (heatmap_width - dis_pairs * (len(genes) - 1)) / (2 * len(genes))
        height_group = heatmap_height / len(groups)

    # determine ratios
    if switch_axis:
        width_ratios = (
            [label_width]
            + [width_group] * 2
            + [dis_pairs, width_group, width_group] * (len(groups) - 1)
            + [
                colorbar_width_spacer,
                colorbar_width,
                size_legend_spacer,
                size_legend_width,
                size_legend_label,
            ]
        )
        height_ratios = [label_height, condition_label] + [height_gene] * len(genes)
    else:
        width_ratios = (
            [label_width]
            + [width_gene] * 2
            + [dis_pairs, width_gene, width_gene] * (len(genes) - 1)
            + [
                colorbar_width_spacer,
                colorbar_width,
                size_legend_spacer,
                size_legend_width,
                size_legend_label,
            ]
        )
        height_ratios = [label_height, condition_label] + [height_group] * len(groups)

    ##### generate the base layout of our plot
    # nrows = 1 for labeling, + number of groups that are to be plotted
    # ncols = 1 for labeling, + number of genes that are to be plotted, + 1 for colorbar, +1 as spacer, +1 as size legend
    fig = figure(figsize=(width, height))
    axs = gridspec.GridSpec(
        nrows=n_rows,
        ncols=n_cols,
        wspace=0.02,
        hspace=0.04,
        width_ratios=width_ratios,
        height_ratios=height_ratios,
    )
    # adjust height of the colorlegend so it does not go the entire length of the figure for a large number of rows
    if n_rows < 8:
        color_legend = fig.add_subplot(axs[2:, -4])
    else:
        colorbar_height = min(5.0, height - label_height)
        wspace = 10.5 / width
        axs2 = gridspec.GridSpecFromSubplotSpec(
            2,
            1,
            subplot_spec=axs[2:, -4],
            wspace=wspace,
            height_ratios=[
                colorbar_height / height,
                (height - colorbar_height) / height,
            ],
        )

        # add ColorBarlegend
        color_legend = fig.add_subplot(axs2[0])

    if plot_absolute:
        max_expression = adata_plot.X.max()
    else:
        max_expression = 1

    cmap = get_cmap("Greys")
    norm = Normalize(vmin=0, vmax=max_expression)
    ColorbarBase(color_legend, cmap=cmap, norm=norm)

    # add legend name (log2cp10k)
    label = fig.add_subplot(axs[1, -4])
    label.axis("off")
    text = label.text(0.5, 0.05, "log2cp10k", ha="center", va="bottom", rotation=0)

    # add dotsize legend relating to cellnumber in group
    label_num_cells = fig.add_subplot(axs[1, -2])
    label_num_cells.axis("off")
    text = label_num_cells.text(
        0.0, 0.05, "number of cells", ha="left", va="bottom", rotation=0
    )

    for n in range(n_bins):
        this_bin = bins.bins.tolist()[n]
        radius = bins.upper_limit_out.get(this_bin)

        bin_plot = fig.add_subplot(axs[n + 2, -2])
        bin_plot.axis("equal")
        sizes = [1]
        bin_plot.pie(
            sizes,
            colors=["white"],
            shadow=False,
            center=(0.5, 0.5),
            radius=radius,
            wedgeprops={"edgecolor": "k", "linewidth": 1},
        )

        bin_label = fig.add_subplot(axs[n + 2, -1])
        bin_label.axis("off")
        text = bin_label.text(
            0.2,
            0.5,
            str(int(bins.upper_limit_in.get(this_bin))),
            ha="left",
            va="center",
            rotation=0,
        )

    # add expression data
    genes = [value for value in genes if value != "NA"]
    groups = [value for value in groups if value != "NA"]

    if switch_axis:
        sub_num_rows = len(genes)
        sub_num_cols = 2 * len(groups) + len(groups) - 1
    else:
        sub_num_rows = len(groups)
        sub_num_cols = 2 * len(genes) + len(genes) - 1

    col_range_cond1 = list(arange(1, sub_num_cols + 1, 3))
    col_range_cond2 = list(arange(2, sub_num_cols + 1, 3))
    colranges = [col_range_cond1, col_range_cond2]

    input_data = {conditions[0]: condition1, conditions[1]: condition2}
    n_cells = {conditions[0]: n_cells_cond1, conditions[1]: n_cells_cond2}

    if switch_axis:
        for c in range(len(conditions)):
            # get datasubset containing only that condition
            adata_subset_plot = input_data[conditions[c]]
            # iterate through the groups (here on the x-axis)
            for i, g in zip(colranges[c], range(len(groups))):
                adata_subset = adata_subset_plot[
                    adata_subset_plot.obs.get(group_by) == groups[g], :
                ].copy()
                # iterate through all of the genes
                for j in range(sub_num_rows):
                    data = _get_expression_table(
                        adata_subset,
                        genes[j],
                        color_map="Greys",
                        max_expression=max_expression,
                    )
                    if data is not None:
                        radius = n_cells[conditions[c]].radius.get(groups[g])
                        ax1 = fig.add_subplot(axs[j + 2, i])
                        ax1.axis("equal")
                        _generate_circle(data, (0.5, 0.5), radius, ax=ax1)
    else:
        for c in range(len(conditions)):
            # get datasubset containing only that condition
            adata_subset_plot = input_data[conditions[c]]
            # iterate through the groups (here on the y-axis)
            for i in range(sub_num_rows):
                adata_subset = adata_subset_plot[
                    adata_subset_plot.obs.get(group_by) == groups[i], :
                ].copy()
                # iterate through all of the genes
                for j, g in zip(colranges[c], range(len(genes))):
                    data = _get_expression_table(
                        adata_subset,
                        genes[g],
                        color_map="Greys",
                        max_expression=max_expression,
                    )
                    if data is not None:
                        radius = n_cells[conditions[c]].radius.get(groups[i])
                        ax1 = fig.add_subplot(axs[i + 2, j])
                        ax1.axis("equal")
                        _generate_circle(data, (0.5, 0.5), radius, ax=ax1)

    # add condition label
    for c in range(len(conditions)):
        for j in colranges[c]:
            label = fig.add_subplot(axs[1, j])
            label.axis("off")
            text = label.text(
                0.5, 0.05, abbrev_labels[c], ha="center", va="bottom", rotation=0
            )

    # add group labels
    if switch_axis:
        for j, g in zip(colranges[0], range(len(groups))):
            label = fig.add_subplot(axs[0, j : (j + 2)])
            label.axis("off")
            text = label.text(
                0.5, 0.05, groups[g], ha="center", va="bottom", rotation=rotation
            )
            label.annotate(
                "",
                xy=(0, 0),
                xycoords="data",
                xytext=(1, 0),
                textcoords="data",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=0."),
            )
    else:
        for j in range(sub_num_rows):
            label = fig.add_subplot(axs[j + 2, 0])
            label.axis("off")
            text = label.text(0.95, 0.5, groups[j], ha="right", va="center")

    # add gene labels
    if switch_axis:
        for j in range(sub_num_rows):
            label = fig.add_subplot(axs[j + 2, 0])
            label.axis("off")
            text = label.text(0.95, 0.5, genes[j], ha="right", va="center")
    else:
        for j, g in zip(colranges[0], range(len(genes))):
            label = fig.add_subplot(axs[0, j : (j + 2)])
            label.axis("off")
            text = label.text(
                0.5, 0.05, genes[g], ha="center", va="bottom", rotation=rotation
            )
            label.annotate(
                "",
                xy=(0, 0),
                xycoords="data",
                xytext=(1, 0),
                textcoords="data",
                arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=0."),
            )
    fig.tight_layout()
    return fig
