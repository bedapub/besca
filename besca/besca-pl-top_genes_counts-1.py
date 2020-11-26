import besca as bc
import matplotlib.pyplot as plt
adata = bc.datasets.pbmc3k_raw()
bc.pl.top_genes_counts(adata)
