from matplotlib import pyplot as plt
import pandas as pd
import importlib
import networkx as nx

def nomenclature_network(tsv):
    """ 
    Based on a tsv config file of a .gmt datasetfile a network is plotted. This network reflects the connection between parent cells and terms.
    Example tsv: 'besca/datasets/genesets/CellNames_scseqCMs6_config.tsv'
    """
    networkx_import = importlib.util.find_spec('networkx')
    pydot_import = importlib.util.find_spec('pydot')

    if (networkx_import or pydot_import) is None:
        raise ImportError(
            "_nomenclature_network.py requires networkx and pydot. Install with pip install pydot and pip install networkx")
        
    # read tsv file 
    df = pd.read_csv(tsv,sep='\t')

    # By default root parents have the entry "None". we need to replace this with its own name so a network per root is created
    roots = df['Parent'] == 'None'
    for row,root in zip(df.iterrows(),roots):
        if root:
            df.at[row[0],'Parent'] = row[1]['Term']
    
    # We create the network with networkx library
    G = nx.from_pandas_edgelist(df, 'Term', 'Parent')
    nx.draw_networkx(G,nx.nx_pydot.pydot_layout(G),font_size=7)
    plt.tight_layout()
    plt.show()
