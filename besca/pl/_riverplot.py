import numpy as np
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns


def riverplot_2categories(adata, categories, palette=None):
    """Generate a riverplot/sanker diagram between two categories.
    parameters
    ----------
    adata: `AnnData`
        anndata object containing data that is to be visualized
    categories: `['str']`
        list of strings identifying the columns that are to be plotted (should be in adata.obs)
    palette: `dict` :
        optional, dict where keys should be keys of the adata.obs[[categories]] and values colors of the node
    returns
    -------
    Figure
        A matplotlib figure element containing the generated plot. To save the figure this plot will need
        to be passed to a parameter and saved in a second step through the fig.savefig() function call.

    Examples
    --------

    >>> # import libraries and dataset
    >>> import besca as bc
    >>> adata = bc.datasets.Baron2016_processed()
    >>> fig = bc.pl.riverplot_2categories(adata,  [ 'assigned_cluster', 'celltype2'])

    .. plot::

        >>> # import libraries and dataset
        >>> import besca as bc
        >>> adata = bc.datasets.Baron2016_processed()
        >>> # define genes
        >>> fig = bc.pl.riverplot_2categories(adata,  [ 'assigned_cluster', 'celltype2'])
        >>> 
    """
    df = adata.obs.copy()
    labels = [x for x in set(df[categories].values.flatten())]

    if palette is None:
        palette = sns.color_palette(
            "deep", len(labels)).as_hex()
    pale = {}
    for idx, x in enumerate(labels):
        pale[x] = palette[idx]
    st_df = None
    i = 0
    st_df = df.groupby([categories[i], categories[i+1]],
                       observed=True).size().reset_index()
    st_df.columns = ['source', 'target', 'count']

    labelSource = df.get(categories[0]).cat.categories
    labelTarget = df[categories[1]].cat.categories

    mappingSource = {}
    for idx, x in enumerate(labelSource):
        mappingSource[x] = idx

    startTarget = len(labelSource)
    mappingTarget = {}
    for idx, x in enumerate(labelTarget):
        mappingTarget[x] = idx + startTarget
    st_df['sourceID'] = st_df['source'].map(mappingSource)
    st_df['targetID'] = st_df['target'].map(mappingTarget)
    data = dict(
        type='sankey', node=dict(
            pad=15, thickness=20, line=dict(color='black', width=0.5),
            label=np.concatenate(
                (labelSource.to_numpy(), labelTarget.to_numpy())),
             color=[pale.get(x) for x in np.concatenate(
                 (labelSource.to_numpy(), labelTarget.to_numpy()))]
        ),
        link=dict(source=st_df['sourceID'],
                  target=st_df['targetID'], value=st_df['count']),
                )
    # creating figure
    fig = go.Figure(dict(data=[data]))
    return(fig)
