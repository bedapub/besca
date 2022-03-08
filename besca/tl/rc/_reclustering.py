from scanpy.preprocessing import filter_genes_dispersion, log1p, regress_out, neighbors
from scanpy.preprocessing import highly_variable_genes as sc_highly_variable_genes
from scanpy.preprocessing import scale as sc_scale
from scanpy.tools import umap, louvain, leiden
from scanpy.tools import pca as sc_pca
from scanpy.plotting import filter_genes_dispersion as plot_filter
from besca._helper import subset_adata as _subset_adata
import sys


def recluster(
    adata,
    celltype,
    celltype_label="leiden",
    min_mean=0.0125,
    max_mean=4,
    min_disp=0.5,
    resolution=1.0,
    regress_out_key=None,
    random_seed=0,
    show_plot_filter=False,
    method="leiden",
    batch_key=None,
    n_shared=2,
):
    """Perform subclustering on specific celltype to identify subclusters.

    Extract all cells that belong to the pre-labeled celltype into a new
    data subset. This datasubset is initialized with the raw data contained in adata.raw. New highly
    variable genes are selected and a new clustering is performed. The function returns the adata
    subset with the new clustering annotation.

    This can be performed on leiden clusters by setting celltype_label = 'leiden' and passing the
    clusters that are to be selected for reclustering as strings or tuple of strings to the parameter
    celltype.

    Parameters
    ----------
    adata:
        the complete AnnData object of the Dataset.
    celltype: `str` or (`str`)
        string identifying the cluster which is to be filtered out, if more than one is to be selected please
        pass them as a tuple not as a list!
    celltype_label: `str` | default = 'leiden'
        string identifying which column in adata.obs will be matching with the celltype argument.
    min_mean: `float` | default = 0.0125
        the minimum gene expression a gene must have to be considered highly variable
    max_mean: `float` | default = 4
        the maximum gene expression a gene can have to be considered highly variable
    min_disp: `float` | default = 0.5
        the minimum dispersion a gene must have to be considered highly variable
    regress_out_key: `list of str` | default = None
        A list of string identifiers of the adata.obs columns that should be regressed out before
        performing clustering. If None then no regress_out is calculated.
    random_seed: `int` | default = 0
        the random seed that is used to produce reproducible PCA, clustering and UMAP results
    show_plot_filter: `bool` | default = False
        boolian value indicating if a plot showing the filtering results for highly variable gene
        detection should be displayed or not
    method: `str` | default = 'leiden'
        clustering method to use for the reclustering of the datasubset. Possible:louvain/leiden
    batch_key: `str` | default = None
        Specify a batch key if the HVG calculation should be done per batch
    n_share: `int` | default = 3
        Divide the nr. of batched by this nr. to get the shared HVGs considered (e.g. >=1/3 of samples)

    Returns
    -------

    AnnData object containing the subcluster annotated with PCA, nearest neighbors, louvain cluster,
    and UMAP coordinates.

    Examples
    --------

    For a more detailed example of the entire reclustering process please refer to the code examples.

    >>> import besca as bc
    >>> import scanpy as sc
    >>> adata = bc.datasets.pbmc3k_processed()
    >>> adata_subset = bc.tl.rc.recluster(adata, celltype=('0', '1', '3', '6'), resolution = 1.3)
    >>> sc.pl.umap(adata_subset, color = ['leiden', 'CD3G', 'CD8A', 'CD4', 'IL7R', 'NKG7', 'GNLY'])

    """
    if not method in ["leiden", "louvain"]:
        raise ValueError("method argument should be leiden or louvain")
    if type(celltype) == str:
        cluster_subset = _subset_adata(adata, adata.obs.get(celltype_label) == celltype)
    elif type(celltype) == tuple:
        filter = adata.obs.get(celltype_label) == "NONE"
        for i in range(len(celltype)):
            filter = filter | (adata.obs.get(celltype_label) == celltype[i])
        cluster_subset = _subset_adata(adata, filter)
    else:
        sys.exit("specify cluster input as a string or tuple")

    cluster_subset.raw = cluster_subset

    # identify highly variable genes
    sc_highly_variable_genes(
        cluster_subset,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        inplace=True,
        batch_key=batch_key,
    )

    if batch_key != None:
        hvglist = cluster_subset.var["highly_variable"].copy()
        hvglist.loc[
            cluster_subset.var["highly_variable_nbatches"]
            >= len(set(cluster_subset.obs[batch_key])) / n_shared,
        ] = True
        cluster_subset.var["highly_variable"] = hvglist.copy()

    if show_plot_filter:
        pl_highly_variable_genes(cluster_subset, show=True)
    print(
        "In total",
        str(sum(cluster_subset.var.highly_variable)),
        "highly variable genes selected within cluster",
    )

    # apply filter
    cluster_subset = _subset_adata(
        cluster_subset, cluster_subset.var.highly_variable, axis=1, raw=False
    )

    # perform further processing
    # log1p(cluster_subset) # data already logged
    if regress_out_key is not None:
        regress_out(cluster_subset, keys=regress_out_key)
    sc_scale(cluster_subset, max_value=10)
    sc_pca(
        cluster_subset, random_state=random_seed, svd_solver="arpack"
    )  # using `svd_solver='arpack' ensures that the PCA leads to reproducible results
    neighbors(cluster_subset, n_neighbors=10, random_state=random_seed)
    umap(cluster_subset, random_state=random_seed)
    if method == "louvain":
        louvain(cluster_subset, resolution=resolution, random_state=random_seed)
    if method == "leiden":
        leiden(cluster_subset, resolution=resolution, random_state=random_seed)

    return cluster_subset


def annotate_new_cellnames(
    adata, cluster_subset, names, new_label="celltype", method="leiden"
):
    """annotate new cellnames to each of the subclusters identified by running recluster.

    Give each subcluster a new celltype identifier. Can only be run on an AnnData subset that has already been
    reclustered, e.g. using the recluster function also available in this package. The list of
    provided names must be of equal length to the number of subclusters in the AnnData subset.
    The order must be the same as the order of the clusters identifiers (integer expected).

    Parameters
    ----------
    adata:
        complete AnnData object from which the subset was extracted. Must contain adata.obs.celltype.
    cluster_subset:
        AnnData subset on which a reclustering has been run. Must contain adata.obs.celltype.
    names: `[str]`
        list of strings containing the new celltype annotation. Order and Length must mach the
        lexicographic order of the louvain clusters in the AnnData susbet.
    celltype_label:
    new_label: `str` | default = 'celltype'
        string specifying under which label in adata.obs the new annotation should be saved
        (it will overwrite existing annotations under this name)
    method: `str` | default = 'leiden'
        string indicating the method used for clustering. This string will be used to retrieve the cluster
        numbers
         in adata.obs and in the cluster_subset. (assuming the same obs column name)

    returns
    -------

    Complete AnnData object with the new celltype annotations in the adata.obs.new_label column for
    the subclusters.

    Examples
    --------

    Please refer to the code examples section of the documentation for a complete example of
    the reclustering process.

    """

    clusters = [
        int(x) for x in cluster_subset.obs.get(method).value_counts().index.tolist()
    ]
    clusters.sort()
    clusters = [str(x) for x in clusters]
    if len(clusters) != len(names):
        print("need to supply as many cluster names as there are clusters")
        print("these are the clusters that need to be annotated:")
        print(str(clusters))
        sys.exit("incorrect number of cluster names supplied")
    else:
        if adata.obs.get(new_label) is None:
            # if the column 'new_label' did not previously exist then create it
            adata.obs[new_label] = "not_labeled"
        else:
            print(
                "NOTE: overwriting labels for the selected cells saved in adata.obs."
                + new_label
                + " with the new labels"
            )
        new_annotation = cluster_subset.obs.get(method).to_frame(name=new_label).copy()

        for i in range(0, len(clusters)):
            new_annotation[new_label].replace(clusters[i], names[i], inplace=True)

        # update anndata object with the new annotation
        adata.obs.update(new_annotation)

        return None
