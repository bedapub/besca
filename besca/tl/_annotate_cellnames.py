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

    # build mapping from old cluster labels to new labels
    label_map = dict(zip(clusters, new_cluster_labels))
    cluster_annotation[new_annotation_label] = cluster_annotation[new_annotation_label].map(label_map)

    # write results back into adata (drop existing column if present to avoid duplicates)
    if new_annotation_label in adata.obs.columns:
        adata.obs = adata.obs.drop(columns=[new_annotation_label])
    adata.obs[new_annotation_label] = cluster_annotation[new_annotation_label]
    return None
