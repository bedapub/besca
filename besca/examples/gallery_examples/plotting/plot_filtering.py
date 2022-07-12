"""
filtering functions
===================

This example shows you how to generate plots to visualize the chosen filter threshold.
This way you can easily check (visually) if your chosen threshold is a good one.

"""

import besca as bc
import matplotlib.pyplot as plt
import pytest
pytest.skip('Test is only for here as example and should not be executed')
adata = bc.datasets.pbmc3k_raw()

# define thresholds
min_genes = 600
min_cells = 2
min_UMI = 600
max_UMI = 6500
max_mito = 0.05
max_genes = 1900

# Visualize filtering thresholds
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(ncols=3, nrows=2)
fig.set_figwidth(15)
fig.set_figheight(8)
fig.tight_layout(pad=4.5)

bc.pl.kp_genes(adata, min_genes=min_genes, ax=ax1)
bc.pl.kp_cells(adata, min_cells=min_cells, ax=ax2)
bc.pl.kp_counts(adata, min_counts=min_UMI, ax=ax3)
bc.pl.max_counts(adata, max_counts=max_UMI, ax=ax4)
bc.pl.max_mito(
    adata, max_mito=max_mito, annotation_type="SYMBOL", species="human", ax=ax5
)
bc.pl.max_genes(adata, max_genes=max_genes)
