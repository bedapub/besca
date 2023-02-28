import besca as bc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import scanpy as sc

def infection_count(adata, col = "n_sgRNAs", experiment="CROPseq"):
    """
    Plots the number of cells against the number of perturbations per cell
    
    parameters
    
    ----------
    adata: `AnnData`
        AnnData object containing data that is to be visualized. Genes need to be in coloumns and cells in rows. Crispr perturbations per cell need to be enumerated.
    col: `str`
        The adata.obs collumn that contains the number of sgRNAs that have perturbed the cell
    experiment: `str`
        Crispr experiment type

    returns
    -------
    None 
        the figure is displayed

    Example
    -------

    Plots the number of cells against the number of perturbations per cell.
    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> adata = bc.crispr_10x_unfiltered()
    >>> bc.pl.infection_count(adata, col = "n_sgRNAs", experiment = "10X")

    .. plot::
        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> import seaborn as sns
        >>> adata = bc.crispr_10x_unfiltered()
        >>> bc.pl.infection_count(adata, col = "n_sgRNAs", experiment = "10X")

    """
    #Get number of cells with all possible sgRNA perturbation counts
    ocur=bc.tl.count_occurrence(adata,col)
    
    #Setup the plot size
    sns.set(font_scale=1.5)
    sns.set(rc={'figure.figsize':(8,5)})
    
    #Barplot with the cell counts on top
    a=sns.barplot(x=ocur.index,y=ocur.Counts/ocur.sum().values,orient='v')
    a.set(xlabel="N. assignments per cell",ylabel="Percentage of total cell count")
    a.set_title(experiment+'\n multiplexity of infection')
    for _,i in ocur.iterrows():
        a.text(i.name, i.Counts/ocur.sum().values, i.Counts, color='black', ha="center")
    
    return

def cell_per_KO(adata, col_target = "Target"):
    """
    Cells per specific KO. If graphs appear to overlap, run again for a better visualization
    
    parameters
    
    ----------
    adata: `AnnData`
        AnnData object containing data that is to be visualized.
    col_target: `str`
        The adata.obs collumn that contains the gene KO per cell
        
    returns
    -------
    None 
        the figure is displayed

    Example
    -------

    Plots Cells per specific KO
    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> adata = bc.crispr_10x_filtered()
    >>> bc.pl.cell_per_KO(adata, col_target = "Target")

    .. plot::
        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> import seaborn as sns
        >>> adata = bc.crispr_10x_filtered()
        >>> bc.pl.cell_per_KO(adata, col_target = "Target")
    """
    # With Control
    with_control=bc.tl.count_occurrence(adata,col_target)
    # Without control
    filtered = filtered[filtered.obs[col_target].str.find("Control") == -1]
    without_control=bc.tl.count_occurrence(filtered,col_target)

    fig, ax =plt.subplots(1,2)
    sns.set(rc={'figure.figsize':(14, adata.obs[col_target].nunique()/3.8)})
    sns.barplot(y=with_control.index,x=with_control.Counts,palette='Blues_d',orient='h', ax=ax[0]).set_title('Number of cells per gene KO (controls included)')
    sns.barplot(y=without_control.index,x=without_control.Counts,palette='Blues_d',orient='h', ax=ax[1]).set_title('Number of cells per gene KO (controls excluded)')
    fig.show()
    plt.subplot_tool()
    
    return

def infection_level(adata, col_target = "Target"):
    """
    Density plot showing infection level of KOs.
    
    parameters
    
    ----------
    adata: `AnnData`
        AnnData object containing data that is to be visualized.
    col_target: `str`
        The adata.obs collumn that contains the gene KO per cell
        
    returns
    -------
    None 
        the figure is displayed
    
    Example
    -------

    Density plot showing infection level of KOs
    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> adata = bc.crispr_10x_filtered()
    >>> bc.pl.infection_level(adata)

    .. plot::
        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> import seaborn as sns
        >>> adata = bc.crispr_10x_filtered()
        >>> bc.pl.infection_level(adata)

    """
    
    occur=bc.tl.count_occurrence(adata,col_target)
    a = sns.displot(data=occur, kde = True, log_scale = True, height=5, aspect=1.5) # make log_scale=True if differences in counts are large
    a.set(xlabel="log(Cell_counts)",ylabel="Gene Target Counts", title='Density cells, mean='+str(occur.mean(axis=0).Counts ))
    
    return None

def _mean_raw(adata):
    if type(adata.X) == np.ndarray:
            adata_X = adata.X
    else:
        adata_X = adata.X.toarray()

    # calculate the percentage of cells that show an expression of the gene above the threshold
    means = np.mean(adata_X, axis=0)

    # add calculated fraction positive to adata.var
    adata.var["mean"] = means.tolist()

def plot_expression_by_sample(adata, gene, col_target ="Target", col_id = "samples__sample_id", raw = True):
    """
    Show mean expression of a gene per sample in a KO
    
    parameters
    
    ----------
    adata: `AnnData`
        AnnData object containing data that is to be visualized.
    col_target: `str`
        The adata.obs collumn that contains the gene KO per cell
    col_id : `str`
        Collumn that contains the sample ids.
    gene: `str` 
        Gene's expression to be ploted
    raw: `bool`
        If the adata we want to use are raw or not
        
    returns
    -------
    None 
        the figure is displayed

    Example
    -------

    Barplot for mean expression of a gene per sample in a KO
    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> adata = bc.crispr_10x_filtered()
    >>> bc.pl.plot_expression_by_sample(adata, "RAB1A", raw = False)

    .. plot::
        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> import seaborn as sns
        >>> adata = bc.crispr_10x_filtered()
        >>> bc.pl.plot_expression_by_sample(adata, "RAB1A", raw = False)
    
    """
    samples = list(adata.obs[col_id].value_counts().index)
    targets = list(adata.obs[col_target].value_counts().index)
    targets.remove('Control')
    per_samples = pd.DataFrame(columns = [col_target, gene + " Expression","Samples"])
    for sample in samples:
        
        filters = adata.obs[col_id] == sample
        temp = bc.subset_adata(adata, filter_criteria= filters, axis = 0, raw=False)
        
        for target in targets:
        
            filters = temp.obs[col_target] == target
            temp_gene = bc.subset_adata(temp, filter_criteria= filters, axis = 0, raw=False)
            if raw:
                    raw_temp = temp_gene.raw
            else:
                raw_temp = temp_gene
            
            # A small variation of the mean function in BESCA that can handle also adata.raw
            _mean_raw(raw_temp)
            expression = raw_temp.var[raw_temp.var.SYMBOL == gene]["mean"][gene]
            # Add the plot only genes that exist and have a non-null expression
            if raw_temp.var[raw_temp.var.SYMBOL == gene] is not None:
                if expression == expression:
                    per_samples = per_samples.append({col_target : target, gene + ' Expression' : expression, 'Samples' : sample}, ignore_index = True)

    plt.figure(figsize=(10, 5))
    sns.boxplot(data = per_samples, x = col_target, y =  gene + " Expression")
    sns.stripplot(data = per_samples, x = col_target, y = gene + " Expression", palette=["black"])
    _ = plt.xticks(rotation=90)

    return None

def avg_KO_persample(adata, col_target = "Target", col_id = "samples__sample_id"):
    """
    Show mean cell count infection of a KO per sample
    
    parameters
    
    ----------
    adata: `AnnData`
        AnnData object containing data that is to be visualized.
    col_target: `str`
        The adata.obs collumn that contains the gene KO per cell
    col_id : `str`
        Collumn that contains the sample ids.
        
    returns
    -------
    None 
        the figure is displayed
    
    Example
    -------
    Dotplot for mean cell count infection of a KO per sample

    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> import numpy as np
    >>> adata = bc.crispr_10x_filtered()
    >>> bc.pl.avg_KO_persample(adata)

    .. plot::
        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> import seaborn as sns
        >>> import numpy as np
        >>> adata = bc.crispr_10x_filtered()
        >>> bc.pl.avg_KO_persample(adata)
    """
    
    counts = pd.DataFrame()
    
    #For each sample get the number of cells perturbed by each guide and merge the rows
    for sample in list(adata.obs[col_id].value_counts().index):
        filters = adata.obs[col_id] == sample
        temp_adata = bc.subset_adata(adata, filter_criteria = filters, raw = False, axis = 0)
        df = bc.tl.count_occurrence(temp_adata, col_target)
        df["sample"] = sample
        counts = counts.append(df)
    
    #Plot the log cell counts
    counts.reset_index(inplace=True)
    counts = counts.rename(columns = {"index" : col_target})
    counts["Counts"] = np.log(counts["Counts"])
    
    _ = plt.xticks(rotation=90)
    sns.scatterplot(data=counts,x = col_target, y = "Counts", hue = "sample")
    plt.title("Perturbations per sample", size=15)
    plt.xlabel("Count of perturbed cells")
    plt.ylabel("Perturbation")
    
    return None

def KO_dotplot(adata, col_target= "Target", title = "Standardized [0,1] expression of targeted genes"):
    """
    Dotplot with expression of each perturbed gene in each perturbation.
    
    parameters
    
    ----------
    adata: `AnnData`
        AnnData object containing data that is to be visualized.
    col_target: `str`
        The adata.obs collumn that contains the gene KO per cell
    title: `str`
        Prefered title for the plot
    returns
    -------
    None 
        the figure is displayed

    Example
    -------
    Dotplot with expression of each perturbed gene in each perturbation.

    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> import scanpy as sc
    >>> adata = bc.crispr_10x_filtered()
    >>> bc.pl.KO_dotplot(adata)
    

    .. plot::
        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> import seaborn as sns
        >>> import scanpy as sc
        >>> adata = bc.crispr_10x_filtered()
        >>> bc.pl.KO_dotplot(adata)
    """
    
    genes=list(adata.obs[col_target][adata.obs[col_target].str.find("Control")==-1].unique())
    genes.sort()
    
    fig_size= adata.obs[col_target].nunique()
    control_last= list(adata.obs[col_target].value_counts().index) #Put the Control catrogory at the end to keep the diagonal of target gene expression and cell with the KO
    control_last.remove('Control')
    control_last.sort()
    control_last.append('Control')
    
    fig = sc.pl.dotplot(adata, var_names= genes , groupby=col_target, title=title, figsize=(2+fig_size/4,1+fig_size/4), show=False, standard_scale = 'var', categories_order = control_last, cmap='viridis' ) 
    fig['mainplot_ax'].set_xlabel('Gene expression')
    fig['mainplot_ax'].set_ylabel('Perturbation identity')

    return

def compute_plot_de_crispr(adata, col_target = "Target", by = None, filter_cells = None):
    """
    Plot the log-fold change of KOs compared to the Control for the QC of Crispr-Cas9 experiments
    
    parameters
    
    ----------
    adata: `AnnData`
        AnnData object containing the sgRNAs per cell.
    col_target: `str`
       The column with the KOs
    by: `str`
        Column that we want to do the dgex with. Example is the samples id columns.
    
    returns
    -------
    None 
        the figure is displayed

    Example
    -------
    Scatterplot of the log-fold change of KOs compared to the Control for the QC of Crispr-Cas9 experiments

    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> adata = bc.crispr_10x_filtered()
    >>> bc.pl.compute_plot_de_crispr(adata)
    

    .. plot::
        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> import seaborn as sns
        >>> adata = bc.crispr_10x_filtered()
        >>> bc.pl.compute_plot_de_crispr(adata)
    """
    # No threshold, so that we don't miss the target genes
    mypval=1
    myfc=0
    dge_data_total = pd.DataFrame()
    if by:
        for sample in list(adata.obs[by].value_counts().index):
            filters = adata.obs[by] == sample
            temp_adata = bc.subset_adata(adata, filter_criteria = filters, raw = False, axis = 0)
            target_list = list(temp_adata.obs[col_target].value_counts().index)
            mypairs = []
            for item in target_list:
                if filter_cells:
                    if (temp_adata.obs[col_target] == item).value_counts()[True]>filter_cells:
                        mypairs.append((item, 'Control'))
                    else:
                        print("Removed " + sample + " for perturbation: " + item +"\n")
                else:
                    mypairs = [(item, 'Control') for item in target_list]
            mypairs.remove(('Control', 'Control'))
            dge_data = bc.tl.crispr.execute_de_sgRNA(temp_adata, mypairs, col_target, myfc, mypval)
            dge_data["sample"] = sample
            dge_data_total = dge_data_total.append(dge_data)
    else:
        conds=list(adata.obs[col_target].cat.categories)
        mypairs = [(item, 'Control') for item in conds]
        mypairs.remove(('Control', 'Control'))
        dge_data_total = bc.tl.crispr.execute_de_sgRNA(adata, mypairs, col_target, myfc, mypval)
    
    cm = plt.cm.get_cmap('viridis')
    
    g = plt.scatter(x=dge_data_total["Name"], y=dge_data_total["Log2FC"],c=dge_data_total["P.adj"], cmap=cm)
    plt.rcParams["figure.figsize"] = (adata.obs[col_target].nunique()/2,5)
    
    plt.colorbar(g, label="-log10(pvalue)")
    _ = plt.xticks(rotation = 90) # Rotates X-Axis Ticks by 90-degrees
    plt.title("Perturbation effects", size=15)
    plt.xlabel("Perturbation")
    plt.ylabel("LogFC control vs KO")
    plt.grid(True)

    return None

def enrichement_per_cluster(enrich_data, col_target = "Target", usr_title = 'KO on each cluster', keep_control = False):
    """
    After Leiden Clustering, plots a heatmap and a batplot of the perturbations counts per cluster
    
    parameters
    
    ----------
    enrich_data: `AnnData`
        AnnData object after leiden clustering.
    col_target: `str`
       The column with the KOs
    usr_tittle: `str`
       Prefered title on the plot
    keep_control: `bool`
        If the control is required
    
    returns
    -------
    None 
        the figure is displayed

    Example
    -------
    Heatmap and a batplot of the perturbations counts per cluster after Leiden

    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> adata = bc.crispr_10x_filtered()
    >>> bc.pl.enrichement_per_cluster(adata, keep_control=True)
    

    .. plot::
        >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> import seaborn as sns
        >>> adata = bc.crispr_10x_filtered()
        >>> bc.pl.enrichement_per_cluster(adata, keep_control=True)
    """
    counts_celltype = bc.tl.count_occurrence_subset(enrich_data, subset_variable = col_target, count_variable='leiden', return_percentage = False)
    if not keep_control:
        counts_celltype = counts_celltype.drop(columns = "Control")

    #Turns counts into percentages
    percentages = pd.DataFrame(index = counts_celltype.index.tolist(), columns = counts_celltype.columns.tolist())
    for cell in counts_celltype.index.tolist():
            data = counts_celltype.loc[cell].tolist()
            percentage = [round(x/sum(data),3) for x in data]
            percentages.loc[cell] = percentage
    
    # Plot as a heatmap and barplot
    fig = percentages.plot(kind='bar', stacked=True, figsize=(10,10))
    fig.set_ylabel('percentage')
    fig.legend(bbox_to_anchor=(1, 1.04), ncol=1, title=usr_title)
    
    fig, ax = plt.subplots(figsize=(11.7,8.27))         # Sample figsize in inches
    fig=sns.heatmap(percentages.astype(float), annot = False, linewidths=.5, ax=ax, cmap="PiYG")
    fig.text(0.5, 1.05,usr_title, fontsize=15,horizontalalignment='center',
     verticalalignment='center',
     transform = ax.transAxes)
    
    return None

def plot_comparison_of_cells(adata, root_path, method="Same", experiment = "CROPseq"):
    """
    Plots the differences between cells of the same or different KO based on the given method. Also stores the data in a file to speed-up execution
    
    parameters
    
    ----------
    adata: `AnnData`
        AnnData object after PCA.
    root_path : `str`
        Path to save the data files
    method: `str`
        Comparison group. Can be "Same", "Different" or "Both"
    experiment: `str`
        Type of Crispr Experiment. Can be CROPseq or 10xChromium
        
    returns
    -------
    A figure is displayed and returns
    the distances of cells based on the method.

    Example
    -------
    Plots the differences between cells of the same or different KO based on the given method

    >>> import besca as bc
    >>> import matplotlib.pyplot as plt
    >>> import seaborn as sns
    >>> import os
    >>> import numpy as np
    >>> adata = bc.crispr_10x_filtered()
    >>> bc.pl.plot_comparison_of_cells(adata, "./", method="Both", experiment = "10xChromium")
    

    .. plot::
         >>> import besca as bc
        >>> import matplotlib.pyplot as plt
        >>> import seaborn as sns
        >>> import os
        >>> import numpy as np
        >>> adata = bc.crispr_10x_filtered()
        >>> bc.pl.plot_comparison_of_cells(adata, "./", method="Both", experiment = "10xChromium")
    """

    #Sanity checks
    if method not in ["Same","Different","Both"]:
        print("Make sure you are giving the correct value for the `method` variable")
        return None

    elif not all(x in adata.obs.columns for x in ['sgRNAs','Target','leiden']):
        print("Make sure adata.obs contains ['sgRNAs','Target','leiden']")
        return None

    if experiment == "CROPseq":
            df=adata.obs[['group', 'sgRNAs','Target','leiden']].copy()
    else:
        df=adata.obs[['sgRNAs','Target','leiden']].copy()

    distances = []
    if method == "Both":
        todos = ["Same", "Different"]
    else:
        todos = [method]
    for todo in todos:
        dist_file = root_path+'sgRNA_distances_'+ todo + '.txt'
        
        #Check if the file exists so that you don't have to do the same again
        if os.path.isfile(dist_file):
            print ("sgRNA distance exists... reading")
            distances.append(pd.read_csv(dist_file, index_col=0, sep="\t"))
        else:
            print("Running sgRNA distance...")
            
            genes=pd.unique(df['sgRNAs'])
            dist_sgRNA=pd.DataFrame(columns=['sgRNAs','Gene','Mean_EucDist'])
            Euc_dist = bc.tl.crispr.find_distances(adata, experiment)
            if todo == "Same":
                for g in genes:
                    # Which rows(cells) are from the same gene 
                    ind=df.index[df['sgRNAs'] == g]
                    # eucledian on pnly the same genes
                    mean_eu=np.nanmean(Euc_dist.loc[ind,ind]) 
                    diag_nr = len(Euc_dist.loc[ind,ind])
                    total_nr = sum(Euc_dist.loc[ind,ind].count())
                    mean_eu_new = (mean_eu * (total_nr - diag_nr) ) / total_nr

                    dist_sgRNA = dist_sgRNA.append({'sgRNAs': g,'Gene': g.split('_')[0], 'Mean_EucDist': mean_eu_new}, ignore_index=True) 
                dist_sgRNA.to_csv(dist_file, sep='\t', header=True)
                distances.append(dist_sgRNA)
            else:
                for g in genes:
                    ind=df.index[df['sgRNAs'] == g]
                    ind_rev=df.index[df['sgRNAs'] != g]
                    # eucledian
                    mean_eu=np.nanmean(Euc_dist.loc[ind,ind_rev]) 

                    dist_sgRNA = dist_sgRNA.append({'gRNA': g, 'Gene': g.split('_')[0], 'Mean_EucDist': mean_eu},ignore_index=True) 
                dist_sgRNA.to_csv(dist_file, sep='\t', header=True)
                distances.append(dist_sgRNA)

    #Plot the boxplots depending on the given method
    fig, ax = plt.subplots(figsize=(30,10)) 
    sns.set(font_scale=1)
    #Sort the means to improve interpretavility
    gene_list = [(distances[0][distances[0]["Gene"] == g].mean(axis=0)["Mean_EucDist"],g) for g in distances[0]["Gene"].value_counts().index]
    gene_list.sort(key = lambda i:i[0], reverse = True)
    indexes = [a for _,a in gene_list]
    if method == "Both":

        distances[0]['type_comparison'] = 'same_gene'
        distances[1]['type_comparison'] = 'diff_gene'
        dist_all = pd.concat(distances)
        sns.boxplot(data=dist_all, x='Gene', y='Mean_EucDist', hue='type_comparison', order = indexes)
    else:

        sns.boxplot(data=distances[0], x='Gene', y='Mean_EucDist', order = indexes)
    ax.set(ylim=(-0.5, max(distances[0]["Mean_EucDist"].max()+2,distances[0]["Mean_EucDist"].max()+2)))
    ax.set_title('Target ' + method + 'type genes')
    ax.set_ylabel('Pairwise EUC distance between cells (PC1,PC2)')
    ax.set_xlabel('Gene')
    return distances