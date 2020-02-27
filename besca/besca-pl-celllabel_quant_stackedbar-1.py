import besca as bc
adata = bc.datasets.pbmc_storage_processed_bbknn()
fig = bc.pl.celllabel_quant_stackedbar(adata, count_variable = 'louvain', subset_variable = 'donor');
