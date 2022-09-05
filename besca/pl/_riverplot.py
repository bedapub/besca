import anndata
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns

def riverplot_2categories(
    adata: anndata.AnnData,
    categories: list,
    palette: dict = None,
    threshold: int = None,
    figsize=None
) -> go:
    """Generate a riverplot/sanker diagram between two categories.
    parameters
    ----------
    adata: `AnnData`
        anndata object containing data that is to be visualized
    categories: `['str']`
        list of strings identifying the columns that are to be plotted (should be in adata.obs)
    palette: `dict` :
        optional, dict where keys should be keys of the adata.obs[[categories]] and values colors of the node
    threshold: `int`
        optional, threshold value, links below threshold will not be display
    figsize: (width, height) or None | default = None
        optional parameter to define the figure size of the plot that is to be generated
    returns
    -------
    Figure
        A matplotlib figure element containing the generated plot. To save the figure this plot will need
        to be passed to a parameter and saved in a second step through the fig.savefig() function call.

    Examples
    --------
    >>> # import libraries and dataset
    >>> import besca as bc
    >>> adata = bc.datasets.simulated_Baron2016_processed()
    >>> fig = bc.pl.riverplot_2categories(adata,  [ 'assigned_cluster', 'celltype2'])

    .. plot::
        >>> # import libraries and dataset
        >>> import besca as bc
        >>> adata = bc.datasets.simulated_Baron2016_processed()
        >>> # define genes
        >>> figure = bc.pl.riverplot_2categories(adata,  [ 'assigned_cluster', 'celltype2'])
        >>> figure = bc.pl.riverplot_2categories(adata,  [ 'assigned_cluster', 'celltype2'], threshold = 15)
    """
    df = adata.obs.copy()
    labels = [x for x in set(df[categories].values.flatten())]

    if palette is None:
        palette = sns.color_palette("deep", len(labels)).as_hex()
    pale = {}
    for idx, x in enumerate(labels):
        pale[x] = palette[idx]
    st_df = None
    i = 0
    st_df = (
        df.groupby([categories[i], categories[i + 1]], observed=True)
        .size()
        .reset_index()
    )
    st_df.columns = ["source", "target", "count"]

    labelSource = df.get(categories[0]).cat.categories
    labelTarget = df[categories[1]].cat.categories

    mappingSource = {}
    for idx, x in enumerate(labelSource):
        mappingSource[x] = idx

    startTarget = len(labelSource)
    mappingTarget = {}
    for idx, x in enumerate(labelTarget):
        mappingTarget[x] = idx + startTarget
    st_df["sourceID"] = st_df["source"].map(mappingSource)
    st_df["targetID"] = st_df["target"].map(mappingTarget)
    # Removing "weak" link for clarity
    if threshold is not None:
        st_df = st_df.loc[st_df.get("count") > threshold]
    data = dict(
        type="sankey",
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=np.concatenate((labelSource.to_numpy(), labelTarget.to_numpy())),
            color=[
                pale.get(x)
                for x in np.concatenate(
                    (labelSource.to_numpy(), labelTarget.to_numpy())
                )
            ],
        ),
        link=dict(
            source=st_df["sourceID"], target=st_df["targetID"], value=st_df["count"]
        ),
    )

    # creating figure
    fig = go.Figure(dict(data=[data]))
    if figsize is not None:
        fig.set_figheight(figsize[1])
        fig.set_figwidth(figsize[0])
    return fig
