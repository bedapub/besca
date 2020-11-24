import besca as bc
import matplotlib.pyplot as plt
adata = bc.datasets.pbmc3k_raw()
min_genes = 600
fig, ax1 = plt.subplots(1)
bc.pl.kp_genes(adata, min_genes = min_genes, ax = ax1)
