"""
Comparing categorical variable
===================

This example shows you how to generate riverplots to compare categorical columns, 
for example to compare multiple annotations
This way you can easily check (visually) discripancies.

"""


import besca as bc

# import data
adata = bc.datasets.Baron2016_processed()

###############################################################################
# compare two categories: annotations made by different annotators
# ----------------------


bc.pl.riverplot_2categories(adata, ["assigned_cluster", "celltype2"])
