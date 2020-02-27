import besca as bc
import matplotlib.pyplot as plt
adata = bc.datasets.brain3k_raw()
bc.pl.librarysize_overview(adata)
