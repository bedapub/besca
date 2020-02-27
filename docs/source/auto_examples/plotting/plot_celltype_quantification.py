"""
visualize cell fractions
========================

This example demonstrates how to generate celltype quantification plots. These types of plots 
can be used to visually represent the number of cells that belong to a certain subset or condition.

"""

import besca as bc 

#import dataset to workwith
adata = bc.datasets.pbmc_storage_processed_bbknn()

#####################
#quantify specific celllabels as a stacked barplot

bc.pl.celllabel_quant_stackedbar(adata, count_variable = 'donor', subset_variable = 'storage_condition');

#####################
# quantify number of cells in each louvain cluster belonging to a specific subset variable (e.g. donor)
#
# Note that the louvain clusters are brought into the correct order.

bc.pl.louvain_quant_stackedbar(adata, subset_variable = 'donor');

#####################
# quantify number of cells belong to each condition in a specific subset
#
# here each dot represents one donor, the boxplots are grouped according to storage condition

bc.pl.celllabel_quant_boxplot(adata, count_variable = 'louvain', subset_variable = 'donor', condition_identifier = 'storage_condition',  plot_percentage = True);

#####################
# here you can also choose to plot total counts instead of percentages

bc.pl.celllabel_quant_boxplot(adata, count_variable = 'louvain', subset_variable = 'donor', condition_identifier = 'storage_condition',  plot_percentage = False);
