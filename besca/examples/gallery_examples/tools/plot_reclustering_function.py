"""
reclustering on specific louvain clusters
=========================================

This example demonstrates who to perform a reclustering on a selected subset of
louvain clusters. You will want to do this for example during the process of celltype
annotation, when the louvain clusters do not have a sufficient resolution to seperate
all clusters and mixed cell populations still exist.

"""

import besca as bc
import scanpy.api as sc

#load and preprocess data (here we will start from a preprocessed dataset)
adata = bc.datasets.pbmc3k_processed()

#extract subset using the recluster function whcih is part of the reclustering (rc) toolkit
adata_subset = bc.tl.rc.recluster(adata, celltype=('0', '1', '3', '6'), celltype_label = 'louvain', resolution = 1.3)

#visualize the new clusters
sc.pl.umap(adata_subset, color = ['louvain', 'CD3G', 'CD8A', 'CD4', 'IL7R', 'NKG7', 'GNLY'])

#append new celltype labels to the subclusters
new_labels = ["CD4 T-cell", #0
              "CD4 T-cell", #1
              "CD4 T-cell", #2
              "CD8 T-cell", #3
              "NK cell", #4
              "CD8 T-cell", #5
              "CD8 T-cell",#6
              "CD4 T-cell", #7
              "CD4 T-cell", #8
              "CD4 T-cell", #9
              "CD4 T-cell"] #10

#merge the labels back into the original adata object
#note this will overwrite what ever was saved in adata.obs.celltype
bc.tl.rc.annotate_new_cellnames(adata, adata_subset, names=new_labels, new_label = 'celltype')

print(adata.obs.celltype.value_counts())
