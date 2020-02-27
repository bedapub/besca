"""
plotting gene expression
========================

This example shows you some of the different plots you can use to plot gene expression.

"""
import besca as bc

#import data
adata = bc.datasets.brain3k_processed()

###############################################################################
# compare two conditions
# ----------------------
#
# You can use the split violin plot to compare gene expression for two different conditions.

bc.pl.gene_expr_split(adata, genes = ['Trem2', 'Nanos1'], split_variable='condition');

###############################################################################
#
# use a stacked split violin plot to compare this for several genes at the same time

#import dataset
adata = bc.datasets.pbmc_storage_processed()

#update the dataset to only contain two conditions
adata = adata[adata.obs.donor != 'Donor_3A']
bc.pl.gene_expr_split_stacked(adata=adata, genes=['CD4', 'CD8A'], split_variable='donor', subset_variable = 'storage_condition');

