import besca as bc
adata = bc.datasets.pbmc_storage_processed_bbknn()
fig = bc.pl.louvain_quant_stackedbar(adata, subset_variable = 'donor');
