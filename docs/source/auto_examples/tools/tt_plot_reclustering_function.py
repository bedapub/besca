"""
reclustering on specific leiden clusters
=========================================

This example demonstrates who to perform a reclustering on a selected subset of
leiden clusters. You will want to do this for example during the process of celltype
annotation, when the leiden clusters do not have a sufficient resolution to seperate
all clusters and mixed cell populations still exist.

"""

import besca as bc
import scanpy as sc

#load and preprocess data (here we will start from a preprocessed dataset)
adata = bc.datasets.pbmc3k_processed()

#extract subset using the recluster function whcih is part of the reclustering (rc) toolkit
adata_subset = bc.tl.rc.recluster(adata, celltype=('CD4-positive, alpha-beta T cell', 'CD8-positive, alpha-beta T cell'), celltype_label = 'celltype2', resolution = 1)

#visualize the new clusters
sc.pl.umap(adata_subset, color = ['leiden',  'CD3G', 'CD8A', 'CD8B','CD4', 'IL7R', 'NKG7', 'GNLY'], color_map = 'viridis')



# We advise to go back to the annotation procedures using auto-annot/sig-annot. 
# As an example here, we performed an a-priori hand annotation.

#append new celltype labels to the subclusters
new_labels = ["CD4 T-cell", #0
              "CD4 T-cell", #1
              "CD8 T-cell", #2
              "NK cell", #3
              "CD8 T-cell", #4
              "CD8 T-cell", #5
              "CD4 T-cell",#6
              "T cell" #7
              ] #10

#merge the labels back into the original adata object
#note this will overwrite what ever was saved in adata.obs.celltype
bc.tl.rc.annotate_new_cellnames(adata, adata_subset, names=new_labels, new_label = 'celltype_rc')


print(adata.obs.celltype_rc.value_counts())
