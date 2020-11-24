# import libraries and dataset
import besca as bc
adata = bc.datasets.Kotliarov2020_processed()
# define genes
genes = ['CD3D', 'CD19']
fig = bc.pl.dot_heatmap(adata, genes=genes, group_by='leiden')
