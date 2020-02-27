import besca as bc
import matplotlib.pyplot as plt
adata = bc.datasets.pbmc3k_raw()
min_cells= 2
fig, ax1 = plt.subplots(1)
bc.pl.kp_cells(adata, min_cells = min_cells, ax = ax1)
