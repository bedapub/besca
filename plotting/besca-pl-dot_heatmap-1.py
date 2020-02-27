#i mport libraries and dataset
import besca as bc
adata = bc.datasets.pbmc_storage_processed()
# define genes
genes = ['CD3D', 'TMSB4X']
fig = bc.pl.dot_heatmap(adata, genes=genes, group_by='louvain')
