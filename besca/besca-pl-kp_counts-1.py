import besca as bc
import matplotlib.pyplot as plt
adata = bc.datasets.pbmc3k_raw()
min_counts = 600
fig, ax1 = plt.subplots(1)
bc.pl.kp_counts(adata, min_counts = min_counts, ax = ax1)
