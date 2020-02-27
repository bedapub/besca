import besca as bc
adata = bc.datasets.pbmc_storage_processed_bbknn()
fig = bc.pl.celllabel_quant_boxplot(adata, count_variable = 'louvain', subset_variable = 'donor', condition_identifier = 'storage_condition',  plot_percentage = True);
