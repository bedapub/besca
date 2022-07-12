import importlib

import networkx as nx
import pandas as pd
from matplotlib import pyplot as plt

def nomenclature_network(
    config_file: str,
    selected_roots=[],
    root_term="None",
    font_size=7,
    node_size=200,
    node_color="tan",
    alpha=0.8,
    figsize=(7,7)
):
    """Plot a nomenclature network based on annotation config file.

    This function plots the relations between celltypes as described within an annotation config file, as the one provided with besca.
    It displays parent - child term relation as a directed graph G ( V, E); Subsetting of such graph is possible using selected_roots argument.


    Parameters
    ----------
    config_file: `str`
        config file from besca, expects a path to a tab separated file containing Parent and Term columns
    selected_roots: `list`
        if list contained terms, will only display the hierarchy starting from those terms.
    root_term : `str`
        the string indicating in the config file that a term does not have a parent term.
    figsize: (width, height) or None | default = (7,7)
        optional parameter to define the figure size of the plot that is to be generated

    Returns
    -------
    Figure
        A matplotlib plt object containing the generated plot.

    Example
    -------
    >>> import besca as bc
    >>> import pkg_resources
    >>> config_file = pkg_resources.resource_filename('besca', 'datasets/genesets/CellNames_scseqCMs6_config.tsv')
    >>> plt = bc.pl.nomenclature_network(config_file)
    >>> plt.show()
    >>> plt = bc.pl.nomenclature_network(config_file, selected_roots = ['Epithelial', 'Tcell'])
    >>> plt.show()

    """
    pydot_import = importlib.util.find_spec("pydot")

    if pydot_import is None:
        raise ImportError(
            "_nomenclature_network.py requires pydot. Install with pip install pydot"
        )
    # read tsv file
    df = pd.read_csv(config_file, sep="\t")

    # By default root parents have the entry "None". we need to replace this with its own name so a network per root is created
    roots_to_set = df["Parent"] == root_term
    for row, root in zip(df.iterrows(), roots_to_set):
        if root:
            df.at[row[0], "Parent"] = row[1]["Term"]

    # We create the network with networkx library
    G = nx.from_pandas_edgelist(
        df, target="Term", source="Parent", create_using=nx.DiGraph()
    )
    ## Subgraph extraction if specific roots were given
    if selected_roots:
        selected_nodes = set()
        for ss in selected_roots:
            try:
                selected_nodes.update(nx.descendants(G, ss))
            except Exception as e:
                print(ss + " node not found in config file")

        G = G.subgraph(list(selected_nodes) + selected_roots)

    if figsize is not None:
        plt.figure(3, figsize=(figsize[0], figsize[1]))
    nx.draw_networkx(
        G,
        nx.nx_pydot.pydot_layout(G),
        font_size=font_size,
        node_size=node_size,
        node_color=node_color,
        alpha=alpha,
    )
    plt.axis("off")
    plt.tight_layout()

    return plt
