#import libraries and dataset
import besca as bc
adata = bc.datasets.pbmc_storage_processed()
#filter out one donor to have two "conditions"
adata = adata[adata.obs.donor != 'Donor_3A',:]
# define genes
genes = ['CD3D', 'TMSB4X']
abbreviations = ['D1', 'D2']
fig = bc.pl.dot_heatmap_split_greyscale(adata, genes=genes, group_by='louvain', split_by = 'donor', , abbreviations = abbreviations)
