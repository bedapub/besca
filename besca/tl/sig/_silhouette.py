import sys
from dataclasses import dataclass

import anndata
import matplotlib
import pandas as pd
import seaborn as sns
import sklearn


@dataclass
class silhouette_in:
    show_samples: matplotlib.axes  # acces is obj.show_samples.get_figure()
    averaged: float


def silhouette_computation(
    adata: anndata.AnnData,
    cluster: str = "dblabel",
    emb: str = "X_umap",
    verbose: bool = False,
) -> silhouette_in:
    """Compute the average and per cell (ie samples) silhouette score for
    the cluster label (should be present in dataobs) (level 3 annotation),
    computed level 2 annotation and a random cell assignbation.
    Return a silhouette_in object

    parameters
    ---------
    adata:  anndata.AnnData
    cluster: 'str'
      clustering to evaluate (should be a column in adata.obs)
    emb: str
        embedding to use for computing the euclidian distance.
        should be a key of  obsm

    returns
    -------
    silhouette_in dataclass object

     Example
    -------
    >>> import besca as bc
    >>> adata = bc.datasets.pbmc3k_processed()
    >>> sils = bc.tl.sig.silhouette_computation (adata)
    >>> print( sils.average)
    >>> sils.show_samples.get_figure()

    """
    if cluster not in adata.obs.keys():
        sys.exit(cluster + " label not found in the dataset (should be in obs)")
    if emb not in adata.obsm.keys():
        sys.exit(emb + " not found in the dataset (should be in obsm)")

    silhouette_avg = sklearn.metrics.silhouette_score(
        adata.obsm[emb], adata.obs.get(cluster)
    )
    sample_silhouette_values = sklearn.metrics.silhouette_samples(
        adata.obsm[emb], adata.obs.get(cluster)
    )
    if verbose:
        print("The average silhouette_score is :", silhouette_avg)
    cluster_labels = adata.obs.get(cluster).unique()
    n_clusters = len(cluster_labels)
    ith_value = {}
    size_cluster = {}
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[
            adata.obs.get(cluster) == cluster_labels[i]
        ]
        ith_value[cluster_labels[i]] = ith_cluster_silhouette_values
        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        size_cluster[cluster_labels[i]] = size_cluster_i

    long_df_silhouette = pd.concat(
        pd.DataFrame({"label": k, "silhouette": v}) for k, v in ith_value.items()
    )
    # TODO ADD PROPER TITLE/ AXIS/ MEAN OVERALL  and mean per label ?
    ax1 = sns.violinplot(
        y=long_df_silhouette["label"], x=long_df_silhouette["silhouette"], scale="count"
    )
    matplotlib.pyplot.close()  # Avoid plooting in function; bad practice
    silhouette_results = silhouette_in(ax1, silhouette_avg)

    return silhouette_results
