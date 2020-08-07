import besca as bc
adata = bc.datasets.Kotliarov2020_processed()
adata.obs   = adata.obs.astype( {'batch' :  'category'})
fig = bc.pl.celllabel_quant_stackedbar(adata, count_variable = 'leiden', subset_variable = 'batch')
