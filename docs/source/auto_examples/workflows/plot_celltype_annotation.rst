.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_workflows_plot_celltype_annotation.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_workflows_plot_celltype_annotation.py:


annotate celltypes
==================

An example workflow using the PBMC3k dataset included with besca illustrating how to annotate celltypes based on Leiden clusters.
This workflow begins with a preprocessed and filtered dataset. 
Please refer to other tutorials on how to perform these steps.
Here we demonstrate a simple annotation based on some known markers.
However, Besca is distributed with auto-annot and sig-annot to annotate automatically or quasi-automatically in a reproducible way datasets.
Please refer to the corresponding tutorials to see how this would work.




.. code-block:: python

    #load libraries
    import besca as bc
    import scanpy as sc

    #load preprocessed dataset (included in BESCA for demonstration purposes)
    adata = bc.datasets.pbmc3k_processed()

    #need to drop celltype annotation stored in this dataset (only relevant for this tutorial)
    adata.obs.drop(columns = ['dblabel', 'leiden'], inplace = True)


    sc.tl.leiden( adata)
    #visualize the louvain clusters
    sc.pl.umap(adata, color=['leiden'])

.. image:: /auto_examples/workflows/images/sphx_glr_plot_celltype_annotation_001.png
    :class: sphx-glr-single-img




visualization of marker genes
-----------------------------

Depending on the type of data you are analysing you will need to look at
different marker genes that are specific to the celltypes you would expect 
to find in your dataset. In this case we are looking at PBMC cells and will
try to identify the main Immunecell subtypes: T-cells, B-cells, Monocytes, and
Dendritic cells.



.. code-block:: python


    #identification of T-cells
    sc.pl.umap(adata, color = ['CD3E', 'CD3G', 'CD3D'])

    #identification of NK cells
    sc.pl.umap(adata, color = ['NCAM1', 'NKG7', 'GNLY'])

    #identification of B-cells
    sc.pl.umap(adata, color = ['MS4A1', 'CD19', 'CD79A'])

    #identification of myeloid cells/dendritic cells
    sc.pl.umap(adata, color = ['CST3', 'S100A8', 'S100A9'])

    #identification of dendritic cells(FCERIA) and monocytes
    sc.pl.umap(adata, color = ['FCER1A','CD14', 'FCGR3A'])




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/workflows/images/sphx_glr_plot_celltype_annotation_002.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/workflows/images/sphx_glr_plot_celltype_annotation_003.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/workflows/images/sphx_glr_plot_celltype_annotation_004.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/workflows/images/sphx_glr_plot_celltype_annotation_005.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/workflows/images/sphx_glr_plot_celltype_annotation_006.png
            :class: sphx-glr-multi-img




cluster level celltype annotation
---------------------------------

Depending on how fine-grained your clustering is you will often come into the
the situation that a leiden cluster contains several cell-populations that
are clearly segregated based on the marker gene expression. If this happens you
can try and adjust the louvain resolution parameter to make the clustering more
fine-grained, but this will not always be successfull. Especially in cases where
your sample contains vastly different celltypes (e.g. an Immuncell dataset 
containing B-cells and T-cells) it will be difficult to resolve T-cell subsets
since they are much more comparable to each other than e.g. a T-cell vs a B-cell.
In this case it often makes sense to make a highlevel cell-labeling and then perform
a second clustering on only the mixed cell clusters. This is the procedure that will
be demonstrated in the rest of this tutorial.



.. code-block:: python


   
    #define high-level celltype annotation
    new_labels = ["T cell", #0
                  "monocyte", #1
                  "mixed", #2
                  "B-cell", #3
                  "T cell", #4
                  "FCGR3A+ monocyte", #5
                  "pDC"] #6

    bc.tl.annotate_cells_clustering(adata, new_labels)

    #visualize annotation
    sc.pl.umap(adata, color = ['celltype'])

    #preserve highlevel labels for future use if desired
    adata.obs['high_level celltype'] = adata.obs.celltype.tolist()






.. image:: /auto_examples/workflows/images/sphx_glr_plot_celltype_annotation_007.png
    :class: sphx-glr-single-img




reclustering on mixed cell clusters
-----------------------------------



.. code-block:: python


    #perform reclustering on subset using besca function 
    adata_subset = bc.tl.rc.recluster(adata, ('mixed', 'T cell'), resolution = 1.2, celltype_label='celltype')
    #visualize important marker genes in reclustering
  
    sc.pl.umap(adata_subset, color = ['leiden', 'CD3G', 'CD8A', 'CD4', 'IL7R', 'NKG7', 'GNLY'], ncols = 3)
 
    
 
    #append new celltype labels to the subclusters
    new_labels = ["CD4 T-cell", #0
                  "CD8 T-cell", #1
                  "CD4 T-cell", #2
                  "CD8 T-cell", #3
                  "NK cell", #4
                  "CD4 T-cell", #5
                  "CD8 T-cell",#6
                  "CD4 T-cell", #7
                  "CD8 T-cell", #8
                  "CD4 T-cell", #9
                  "CD4 T-cell", #10
                  "CD4 T-cell" ] #11
    
    #merge the labels back into the original adata object
    #note this will overwrite what ever was saved in adata.obs.celltype
    bc.tl.rc.annotate_new_cellnames(adata, adata_subset, names=new_labels, new_label = 'celltype')
    
    sc.pl.umap(adata, color = ['celltype'])
    print(adata.obs.celltype.value_counts())


.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/workflows/images/sphx_glr_plot_celltype_annotation_008.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/workflows/images/sphx_glr_plot_celltype_annotation_009.png
            :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    In total 1471 highly variable genes selected within cluster
    NOTE: overwriting labels for the selected cells saved in adata.obs.celltype with the new labels


**Total running time of the script:** ( 0 minutes  15.606 seconds)


.. _sphx_glr_download_auto_examples_workflows_plot_celltype_annotation.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: plot_celltype_annotation.py <plot_celltype_annotation.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: plot_celltype_annotation.ipynb <plot_celltype_annotation.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.readthedocs.io>`_
