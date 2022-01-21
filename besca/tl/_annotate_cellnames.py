import sys


def annotate_cells_clustering(
    adata,
    new_cluster_labels,
    new_annotation_label="celltype",
    clustering_label="leiden",
):
    """Function to add annotation to adata.obs based on clustering
    This function replaces the original cluster labels located in the column clustering_label with
    the new values specified in the list new_cluster_lables. The values in this list need to be in the
    same lexicographic order as the cluster values.
    The order of this can be checked by running
    ` sorted(adata.obs.get(clustering_label).unique(), key=int)` or
    ` sorted(adata.obs.get(clustering_label).unique(), key=str)`

    parameters
    ----------

    adata: AnnData
        the AnnData object that is supposed to recieve a new annotation
    new_cluster_labels: `list`
        a list in correct lexicographic order of the new cluster labels
    new_annotation_label: `str` | default = 'celltype'
        string identifying the name underwhich the new cluster labels should be added to adata.obs
    clustering_label: `str` | default = 'louvain'
        string identifying the name underwhich the old cluster labels can be found in adata.obs

    returns
    -------

    AnnData
        The adata.obs has been updated with a new column containing the new annotation

    """
    # get the cluster ids
    if clustering_label not in adata.obs.columns:
        sys.exit(clustering_label + " label not found")
    try:
        clusters = sorted(adata.obs.get(clustering_label).unique(), key=int)
    # If int conversion does not work
    except:
        clusters = sorted(adata.obs.get(clustering_label).unique(), key=str)
    # extract old clustering labels to a new dataframe
    cluster_annotation = (
        adata.obs.get(clustering_label).to_frame(name=new_annotation_label).copy()
    )
    # check if the number of new labels matches the number of clusters
    if len(clusters) != len(new_cluster_labels):
        sys.exit(
            "Specified "
            + str(len(new_cluster_labels))
            + " new labels for a total of "
            + str(len(clusters))
            + " clusters. Numbers should match! No changes were made.\n"
        )

    # replace old labels with the new labels
    for i in range(0, len(clusters)):
        cluster_annotation.get(new_annotation_label).replace(
            clusters[i], new_cluster_labels[i], inplace=True
        )

    # write results back into adata
    adata.obs = adata.obs.merge(
        cluster_annotation, how="outer", left_index=True, right_index=True
    )
    return None
