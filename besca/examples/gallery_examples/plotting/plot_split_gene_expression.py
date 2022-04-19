"""
# TODO
plotting gene expression
========================

This example shows you some of the different plots you can use to plot gene expression.

"""
import besca as bc
import pytest
pytest.skip('Test is only for here as example and should not be executed')
# import data
adata = bc.datasets.Haber2017_processed()

###############################################################################
# compare two conditions
# ----------------------
#
# You can use the split violin plot to compare gene expression for two different conditions.

bc.pl.gene_expr_split(adata, genes=["Defa24", "Gm15284"], split_variable="donor")

###############################################################################
#
# use a stacked split violin plot to compare this for several genes at the same time

bc.pl.gene_expr_split_stacked( # ISSUE is in this method
    adata=adata,
    genes=["Defa24", "Gm15284"],
    split_variable="donor",
    subset_variable="region_x",
)
