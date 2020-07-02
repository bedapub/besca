from matplotlib import pyplot as plt
import pandas as pd
import importlib
import networkx as nx

def nomenclature_network(config_file,mywidth=7,myheight=7,font_size=7,
                        node_size=200,node_color='tan',alpha=0.8,root_term='None'):
    """Plot a nomenclature network based on annotation config file.

    This function plots a nomeclature network based on the annotation config file provided with besca.

    Parameters
    ----------
    config_file: `str`

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

    """
    pydot_import = importlib.util.find_spec('pydot')

    if pydot_import is None:
        raise ImportError(
            "_nomenclature_network.py requires pydot. Install with pip install pydot")
        
    # read tsv file 
    df = pd.read_csv(config_file,sep='\t')

    # By default root parents have the entry "None". we need to replace this with its own name so a network per root is created
    roots = df['Parent'] == root_term
    for row,root in zip(df.iterrows(),roots):
        if root:
            df.at[row[0],'Parent'] = row[1]['Term']
    
    # We create the network with networkx library
    G = nx.from_pandas_edgelist(df, 'Term', 'Parent')
    plt.figure(3,figsize=(mywidth,myheight)) 
    nx.draw_networkx(G,nx.nx_pydot.pydot_layout(G),font_size=font_size,node_size=node_size,node_color=node_color,alpha=alpha)
    plt.axis('off')
    plt.tight_layout()

    return plt
