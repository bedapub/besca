# import libraries and dataset
import besca as bc
adata = bc.datasets.Baron2016_processed()
# define genes
fig = bc.pl.riverplot_2categories(adata,  [ 'assigned_cluster', 'celltype2'])
# >>>
