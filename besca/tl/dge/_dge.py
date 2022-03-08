# general libraries
import os
from scanpy.tools import rank_genes_groups
from pandas import DataFrame, read_csv, concat
from numpy import sign, log10, inf, max, abs, floor, arange, round, sort
from plotly.offline import plot as plotly_plot

# other besca functions
from besca._helper import get_raw, subset_adata


def extract_info_rank_genes_groups(adata):
    """This function extracts information from adata.uns[rank_gene_groups] to pandas.DataFrames

    This function outputs 4 pandas.DataFrames which contain the scores, pvalues, logFC, and FDRs
    calculated by scanpy.tl.rank_gene_groups() and saved to adata.uns['rank_gene_groups'].

    This function can only be run on AnnData objects that have already been passed to scanpy.tl.rank_gene_groups().

    parameters
    ----------
    adata: `AnnData`
        AnnData object which contains information in adata.uns['rank_gene_groups']

    returns
    -------
    (scores, pvalues, logFC, FDRs)
        each of type pandas.DataFrame

    """
    rank_genes = adata.uns["rank_genes_groups"]
    groups = rank_genes["names"].dtype.names
    group_number = len(groups)

    # get gene information
    mydat = adata.var.loc[:, ["SYMBOL", "ENSEMBL"]]
    mydat.rename(columns={"SYMBOL": "Description"}, inplace=True)

    # initialize dataframe for storing the scores
    scores = DataFrame(
        {"n_" + key[:1]: rank_genes[key][groups[0]] for key in ["names", "scores"]}
    )
    scores.rename(columns={"n_s": groups[0], "n_n": "NAME"}, inplace=True)
    scores.set_index("NAME", inplace=True)
    scores.index = scores.index.astype(str)

    for i in range(1, group_number):
        n_scores = DataFrame(
            {"n_" + key[:1]: rank_genes[key][groups[i]] for key in ["names", "scores"]}
        )
        n_scores.set_index("n_n", inplace=True)
        n_scores.index = n_scores.index.astype(str)
        scores = scores.merge(n_scores, how="left", left_index=True, right_index=True)
        newname = groups[i]
        scores.rename(columns={"n_s": newname}, inplace=True)
    scores = scores.astype(float)

    # merge in gene annotation
    scores = mydat.merge(scores, how="right", left_index=True, right_index=True)

    # make index into ENSEMBL instead of symbol
    scores.set_index("ENSEMBL", inplace=True)
    scores.index.names = ["NAME"]

    # get pvalues
    # initialize dataframe for storing the pvalues
    pvalues = DataFrame(
        {"n_" + key[:1]: rank_genes[key][groups[0]] for key in ["names", "pvals"]}
    )
    pvalues.rename(columns={"n_p": groups[0], "n_n": "NAME"}, inplace=True)
    pvalues.set_index("NAME", inplace=True)
    pvalues.index = pvalues.index.astype(str)

    for i in range(1, group_number):
        n_pvalues = DataFrame(
            {"n_" + key[:1]: rank_genes[key][groups[i]] for key in ["names", "pvals"]}
        )
        n_pvalues.set_index("n_n", inplace=True)
        n_pvalues.index = n_pvalues.index.astype(str)
        pvalues = pvalues.merge(
            n_pvalues, how="left", left_index=True, right_index=True
        )
        newname = groups[i]
        pvalues.rename(columns={"n_p": newname}, inplace=True)
    pvalues = pvalues.astype(float)

    # merge in pvalues
    pvalues = mydat.merge(pvalues, how="right", left_index=True, right_index=True)

    # make index into ENSEMBL instead of symbol
    pvalues.set_index("ENSEMBL", inplace=True)
    pvalues.index.names = ["NAME"]

    # get logFC
    logFC = DataFrame(
        {
            "n_" + key[:1]: rank_genes[key][groups[0]]
            for key in ["names", "logfoldchanges"]
        }
    )
    logFC.rename(columns={"n_l": groups[0], "n_n": "NAME"}, inplace=True)
    logFC.set_index("NAME", inplace=True)
    logFC.index = logFC.index.astype(str)

    for i in range(1, group_number):
        n_logFC = DataFrame(
            {
                "n_" + key[:1]: rank_genes[key][groups[i]]
                for key in ["names", "logfoldchanges"]
            }
        )
        n_logFC.set_index("n_n", inplace=True)
        n_logFC.index = n_logFC.index.astype(str)
        logFC = logFC.merge(n_logFC, how="left", left_index=True, right_index=True)
        newname = groups[i]
        logFC.rename(columns={"n_l": newname}, inplace=True)
    logFC = logFC.astype(float)

    # merge in pvalues
    logFC = mydat.merge(logFC, how="right", left_index=True, right_index=True)

    # make index into ENSEMBL instead of symbol
    logFC.set_index("ENSEMBL", inplace=True)
    logFC.index.names = ["NAME"]

    # get FDRs
    FDRs = DataFrame(
        {"n_" + key[:1]: rank_genes[key][groups[0]] for key in ["names", "pvals_adj"]}
    )
    FDRs.rename(columns={"n_p": groups[0], "n_n": "NAME"}, inplace=True)
    FDRs.set_index("NAME", inplace=True)
    FDRs.index = FDRs.index.astype(str)

    for i in range(1, group_number):
        n_FDRs = DataFrame(
            {
                "n_" + key[:1]: rank_genes[key][groups[i]]
                for key in ["names", "pvals_adj"]
            }
        )
        n_FDRs.set_index("n_n", inplace=True)
        n_FDRs.index = n_FDRs.index.astype(str)
        FDRs = FDRs.merge(n_FDRs, how="left", left_index=True, right_index=True)
        newname = groups[i]
        FDRs.rename(columns={"n_p": newname}, inplace=True)
    FDRs = FDRs.astype(float)

    # merge in pvalues
    FDRs = mydat.merge(FDRs, how="right", left_index=True, right_index=True)

    # make index into ENSEMBL instead of symbol
    FDRs.set_index("ENSEMBL", inplace=True)
    FDRs.index.names = ["NAME"]

    return (scores, pvalues, logFC, FDRs)


def perform_dge(
    adata,
    design_matrix,
    differentiating_criteria,
    constant_criteria,
    basepath,
    min_cells_per_group=30,
    method="wilcoxon",
):
    """Perform differential gene expression between two conditions over many adata subsets.

    This function automatically generates top_tables and rank_files for a list of comparisons
    in a dataset. The comparison you wish to perform need to be identified in a so called design
    matrix (see below).

    This function is capable of handling comparisons where you wish to compare
    two conditions in a subset of the dataset, e.g. treatment vs control in the celltype CD4 T-cell.
    The conditions must be annotated in a column adata.obs, this represents the differentiating_criteria.
    This column may only have two different labels! The subsets in which this comparison should be made
    must be annotated in another column represented by 'constant_criteria'. This column may have as
    many labels as you wish.

    Design Matrix:
    the design matrix consists of a pandas.DataFrame with two columns. Each row
    represents one comparison that is to be made. The first column labeled 'Group1',
    contains a tuple identifying the first group for that comparison and the
    second column labeled 'Group2' contains a tuple identifying the second group for
    the comparison. The tuple has the form (differentiating_criteria, constant_criteria).

    >>> #example of a Design Matrix
    >>> celltypes = ['CD4 T-cell', 'CD8 T-cell', 'B-cell', 'myeloid cell']
    >>> design_matrix = pd.DataFrame({'Group1':[('PBMC', celltype) for celltype in celltypes], 'Group2':[('Skin', celltype) for celltype in celltypes]})

    parameters
    ----------
    adata: `AnnData`
        AnnData object containing
    design_matrix: `pandas.DataFrame`
        pandas.DataFrame containing all the comparisons that are to be made.
    method: `str`
        one of 't-test', 'wilcoxon', 't-test_overestim_var', 'logreg'
    """

    # get raw data
    adata_raw = get_raw(adata=adata)

    too_few_cells = []
    for i in range(design_matrix.shape[0]):

        # get group1 and group2
        group1 = design_matrix.values[i, 0]
        group2 = design_matrix.values[i, 1]

        contrast = str(
            group1[1] + "_" + group1[0] + "_vs_" + group2[1] + "_" + group2[0]
        )
        contrast = contrast.replace(" ", "_")
        contrast = contrast.replace("-", "_")
        contrast = contrast.replace("+", "")

        # perform sanitychecks
        if group1[1] != group2[1]:
            # sys.exit()
            print("please ensure design matrix is correct!")

        # get adata_subset
        adata_subset = subset_adata(
            adata_raw,
            filter_criteria=adata.obs.get(constant_criteria) == group1[1],
            raw=False,
        )
        adata_subset.var = adata_raw.var

        # ensure that you have enough cells in each group
        counts = adata_subset.obs.get(differentiating_criteria).value_counts()
        if counts.get(group1[0]) is None:
            print("----------------------------------------")
            print("not enough cells for comparison", contrast, " in one of the groups")
            too_few_cells.append(i)
        elif counts.get(group2[0]) is None:
            print("----------------------------------------")
            print("not enough cells for comparison", contrast, " in one of the groups")
            too_few_cells.append(i)
        elif not (
            (counts.get(group1[0]) > min_cells_per_group)
            & (counts.get(group2[0]) > min_cells_per_group)
        ):
            print("----------------------------------------")
            print("not enough cells for comparison", contrast, " in one of the groups")
            too_few_cells.append(i)
        else:
            # caclulate rank_genes
            print("----------------------------------------")
            print("performing comparison", contrast)
            rank_genes_groups(
                adata=adata_subset,
                groupby=differentiating_criteria,
                reference=group2[0],
                groups=[group1[0]],
                method=method,
                n_genes=adata_subset.n_vars,
            )

            scores, pvalues, logFC, FDRs = extract_info_rank_genes_groups(adata_subset)

            # calculate mean expression of each gene in the data
            mean_all = adata_subset.X.todense().mean(axis=0).tolist()[0]
            mean_group1 = (
                adata_subset[
                    adata_subset.obs.get(differentiating_criteria) == group1[0]
                ]
                .X.todense()
                .mean(axis=0)
                .tolist()[0]
            )
            mean_group2 = (
                adata_subset[
                    adata_subset.obs.get(differentiating_criteria) == group2[0]
                ]
                .X.todense()
                .mean(axis=0)
                .tolist()[0]
            )

            # get mean expression
            log2cp10k = DataFrame(
                data={
                    "log2cp10k": mean_all,
                    "log2cp10k_" + group1[0]: mean_group1,
                    "log2cp10k_" + group2[0]: mean_group2,
                },
                index=adata_subset.var.ENSEMBL,
            )

            # generate top_table for group1
            top_table = DataFrame(
                data={
                    "Contrast": contrast,
                    "SYMBOL": scores.Description,
                    "ENSEMBL": scores.index.tolist(),
                    "Score": scores.get(group1[0]),
                },
                index=scores.index,
            )
            # add p-value, logFC, FDRs, mean_expression
            top_table = top_table.merge(
                pvalues.get(group1[0])
                .to_frame()
                .rename(columns={group1[0]: "P-value"}),
                how="left",
                right_index=True,
                left_index=True,
            )
            top_table = top_table.merge(
                logFC.get(group1[0]).to_frame().rename(columns={group1[0]: "logFC"}),
                how="left",
                right_index=True,
                left_index=True,
            )
            top_table = top_table.merge(
                FDRs.get(group1[0]).to_frame().rename(columns={group1[0]: "FDR"}),
                how="left",
                right_index=True,
                left_index=True,
            )
            top_table = top_table.merge(
                log2cp10k, how="left", right_index=True, left_index=True
            )

            # order according to score
            top_table.sort_values(ascending=False, axis=0, by="Score", inplace=True)

            # replace all 0 values with the smallest number possible in python
            smallest_num = 1e-308
            top_table.replace(0, 1e-308, inplace=True)

            # generate rank files
            rank_file = top_table.get(["SYMBOL", "P-value", "logFC"])
            rank_file["value"] = abs(log10(rank_file["P-value"])) * sign(
                rank_file["logFC"]
            )
            rank_file.drop(columns=["P-value", "logFC"], inplace=True)
            rank_file.sort_values(ascending=False, axis=0, by="value", inplace=True)

            # replace all inf values with very small or very large values
            rank_file.replace(inf, 1e100, inplace=True)
            rank_file.replace(-inf, -1e100, inplace=True)

            outpath = os.path.join(basepath)
            if not os.path.exists(outpath):
                os.makedirs(outpath)

            if method == "wilcoxon":
                rank_File = os.path.join(outpath, "WilxRank_" + contrast + ".rnk")
                TopTable_File = os.path.join(
                    outpath, "WilxTopTable_" + contrast + ".txt"
                )
            elif method == "t-test_overestim_var":
                rank_File = os.path.join(outpath, "OverestRank_" + contrast + ".rnk")
                TopTable_File = os.path.join(
                    outpath, "OverestTopTable_" + contrast + ".txt"
                )
            elif method == "t-test":
                rank_File = os.path.join(outpath, "tTestRank_" + contrast + ".rnk")
                TopTable_File = os.path.join(
                    outpath, "tTestTopTable_" + contrast + ".txt"
                )
            elif method == "logreg":
                rank_File = os.path.join(outpath, "logregRank_" + contrast + ".rnk")
                TopTable_File = os.path.join(
                    outpath, "logregTopTable_" + contrast + ".txt"
                )
            else:
                sys.exit(
                    "need to specify type as one of 'wilcoxon' or 't-test_overestim_var'  or 't-test' or 'logreg'"
                )

            # write out rankfile
            rank_file.to_csv(
                rank_File, sep="\t", index=False, header=False, float_format="%.2f"
            )
            print(rank_File, "written out")

            # write out TopTable
            top_table.to_csv(
                TopTable_File, sep="\t", index=True, header=True, float_format="%1.3e"
            )
            print(TopTable_File, "written out")

    # return table containing the comparisions that could not be performed
    if len(too_few_cells) > 0:
        print("")
        print("")
        print("----------------------------------------------")
        print(
            "the following comparisons from the design_matrix could not per performed due to too few cells"
        )
        display(design_matrix.iloc[too_few_cells, :])

    return None
    sys.exit(0)


def plot_interactive_volcano(top_table_path, outdir):
    """plot an interactive volcano plot based on toptable file.

    This function takes as input a toptable file generated with the functions x, y
    and uses the contained information to generate an interactive html plot.

    parameters
    ----------
    top_table_path: `str`
        string indicating the filepath to the toptable file
    outdir: `str`
        string indicating the filepath where the generated file is to be written to.

    returns
    -------
    None
        writes an interactive html file out
    """
    # generate proper filename and plot title from filepath
    comparision_name = os.path.basename(top_table_path)
    comparision_name = comparision_name.split("TopTable_", 1)[1]
    comparision_name = comparision_name.strip(".txt")

    filename = "volcano_plot_" + comparision_name + ".html"
    title = comparision_name.replace("_", " ")

    # read in data
    df = read_csv(top_table_path, sep="\t")

    # get group name of expression
    group1 = df.columns[9]
    group1_label = group1.replace("log2cp10k_", "")
    group2 = df.columns[10]
    group2_label = group2.replace("log2cp10k_", "")

    def generate_data_trace(df, threshold):
        """Helper function to iteratively generate the different traces that are plotted."""
        df_filtered = df[df.log2cp10k > threshold]

        # get annotation per point
        gene = ["Gene Symbol: " + gene + "<br>" for gene in df_filtered.SYMBOL.tolist()]
        expression = [
            "log2cp10k: " + str(expr) + "<br>"
            for expr in df_filtered.log2cp10k.tolist()
        ]
        expression_group1 = [
            "log2cp10k of " + group1_label + ": " + str(expr) + "<br>"
            for expr in df_filtered.get(group1).tolist()
        ]
        expression_group2 = [
            "log2cp10k of " + group2_label + ": " + str(expr) + "<br>"
            for expr in df_filtered.get(group2).tolist()
        ]
        pvalue = ["P-value: " + str(expr) + "<br>" for expr in df_filtered.FDR.tolist()]
        foldchange = [
            "logFC: " + str(value) + "<br>" for value in df_filtered.logFC.tolist()
        ]
        text = [
            a + b + c + d + e + f
            for a, b, c, d, e, f in zip(
                gene,
                expression,
                expression_group1,
                expression_group2,
                pvalue,
                foldchange,
            )
        ]

        # generate list of transforms
        minlogFC = min(df_filtered.logFC)
        maxlogFC = max(df_filtered.logFC)

        transforms = [
            dict(type="filter", target="y", operation="<", value=1),
            dict(type="filter", target="x", operation=">", value=minlogFC),
        ]

        # generate data
        data = dict(
            type="scatter",
            x=abs(df_filtered.logFC.tolist()),
            y=df_filtered.FDR.tolist(),
            text=text,
            hoverinfo="text",
            mode="markers",
            opacity=0.9,
            visible=threshold == 0,
            marker=dict(
                size=4,
                color=df_filtered.log2cp10k.tolist(),  # set color equal to a variable
                colorscale="Viridis",
                colorbar=dict(title="log2cp10k"),
                showscale=True,
            ),
            transforms=transforms,
        )
        return data

    # define range of pvalues to generate traces for
    FDRs = [1, 0.05, 0.01, 1e-5, 1e-10, 1e-50, 1e-100]

    # define cutoffpoints for log2cp10k values to display (dynamically depending on range in dataset)
    max_log2cp10k = int(floor(max(df.log2cp10k)))
    log2cp10k = list(round(arange(0.0, 1.5, 0.100000000000), 2)) + list(
        range(2, max_log2cp10k + 1, 1)
    )

    # define cutoffpoints for logFC
    logFC = [0, 1, 2, 3, 4, 5]

    # generate individual trace for each of the log2cp10k thresholds
    traces = []
    for i in log2cp10k:
        traces.append(generate_data_trace(df, i))

    # helper functions to generate the buttons

    def generate_pvalue_button(value):
        button = {
            "method": "restyle",
            "args": ["transforms[0].value", value],
            "label": str(value),
        }
        return button

    def generate_logFC_button(value):
        button = {
            "method": "restyle",
            "args": ["transforms[1].value", value],
            "label": str(value),
        }
        return button

    def generate_log2cp10k_button(value, number, total):
        visible = [False] * total
        visible[number] = True
        button = {
            "method": "restyle",
            "args": [{"visible": visible}],
            "label": str(value),
        }
        return button

    # generate buttons
    pvalue_buttons = []
    for p in FDRs:
        pvalue_buttons.append(generate_pvalue_button(p))

    logFC_buttons = []
    for l in logFC:
        logFC_buttons.append(generate_logFC_button(l))

    cp10k_buttons = []
    for c in range(len(log2cp10k)):
        cp10k_buttons.append(generate_log2cp10k_button(log2cp10k[c], c, len(log2cp10k)))

    updatemenus = list(
        [
            dict(
                buttons=list(pvalue_buttons),
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.2,
                xanchor="left",
                y=1.07,
                yanchor="top",
            ),
            dict(
                buttons=list(logFC_buttons),
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.55,
                xanchor="middle",
                y=1.07,
                yanchor="top",
            ),
            dict(
                buttons=list(cp10k_buttons),
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.85,
                xanchor="middle",
                y=1.07,
                yanchor="top",
            ),
        ]
    )

    ymax = min(df.FDR)
    xmax = max(abs(df.logFC))

    layout = dict(
        xaxis=dict(
            autorange=True,
            range=[0, xmax],
            title="abs(logFC)",
        ),
        yaxis=dict(
            type="log",
            autorange="reversed",
            range=[1, 0],
            title="FDR",
        ),
        title=title,
        hovermode="closest",
        updatemenus=updatemenus,
    )

    annotations = list(
        [
            dict(
                text="FDR: <br> max value",
                x=0.1,
                y=1.05,
                align="left",
                yref="paper",
                xref="paper",
                showarrow=False,
            ),
            dict(
                text="absolute logFC: <br> min value",
                x=0.4,
                y=1.05,
                align="left",
                yref="paper",
                xref="paper",
                showarrow=False,
            ),
            dict(
                text="log2cp10k expression: <br> min value",
                x=0.76,
                y=1.05,
                align="left",
                yref="paper",
                xref="paper",
                showarrow=False,
            ),
        ]
    )

    layout["annotations"] = annotations

    # ensure that basepath folder exists
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print("plotting", filename)
    plotly_plot(
        {"data": traces, "layout": layout},
        validate=False,
        filename=os.path.join(outdir, filename),
        auto_open=False,
    )


def get_de(adata, mygroup, demethod="wilcoxon", topnr=5000, logfc=1, padj=0.05):
    """Get a table of significant DE genes at certain cutoffs
    Based on an AnnData object and an annotation category (e.g. louvain) runs
    scanpy's rank_genes_groups using a specified method with specified cutoffs
    (nr. genes, logfc, padj) and returns a df with the results
    parameters
    ----------
    adata: `AnnData`
        AnnData object containing
    mygroup: `str`
        group for performing de, needs be in adata.obs
    demethod: `str`
        one of 't-test', 'wilcoxon', 't-test_overestim_var', 'logreg'
    topnr: `int`
        the number of top genes in the DE analysis
    padj: `float`
        log fold-change cutoff
    logfc: `float`
        adjusted p-value cutoff
    returns
    -------
    delist
        a list of panda DataFrames of differentially expressed genes
    """

    try:
        x = adata.obs[mygroup]
    except KeyError:
        print(
            "Oops!  The adata object does not have the specified column. Options are: "
        )
        print(list(adata.obs.columns))
        return

    mygroups = list(sort(list(set(adata.obs[mygroup]))))
    delist = {}
    rank_genes_groups(
        adata,
        groupby=mygroup,
        use_raw=True,
        n_genes=adata.raw.X.shape[1],
        method=demethod,
    )
    for i in mygroups:
        df = DataFrame(adata.uns["rank_genes_groups"]["names"]).head(topnr)[i]
        dfS = DataFrame(adata.uns["rank_genes_groups"]["scores"]).head(topnr)[i]
        dfFC = DataFrame(adata.uns["rank_genes_groups"]["logfoldchanges"]).head(topnr)[
            i
        ]
        dfp = DataFrame(adata.uns["rank_genes_groups"]["pvals_adj"]).head(topnr)[i]
        d = concat([df, dfS, dfFC, dfp], axis=1)
        d.columns = ["Name", "Score", "Log2FC", "P.adj"]
        delist[i] = d[(d["Log2FC"] >= logfc) & (d["P.adj"] <= padj)]
    return delist
