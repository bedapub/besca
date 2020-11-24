# import libraries and dataset
import besca as bc
adata = bc.datasets.Haber2017_processed()
# define genes
genes = ['Defa22', 'Defa24', 'Gm15284', 'Reg4']
fig = bc.pl.dot_heatmap_split(adata, genes=genes, group_by='leiden', split_by = 'donor')
