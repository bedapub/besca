import besca as bc
import matplotlib.pyplot as plt
adata = bc.datasets.pbmc3k_raw()
max_counts = 6500
fig, ax1 = plt.subplots(1)
bc.pl.max_counts(adata, max_counts = max_counts, ax = ax1)
