.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_preprocessing_plot_pca_neighbors_clustering.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_auto_examples_preprocessing_plot_pca_neighbors_clustering.py:


cluster generation
==================

This example demonstrates how to perform highly variable gene selection, PCA, nearest neighbor calculation, and clustering.


.. code-block:: default


    import besca as bc
    import scanpy as sc

    #import example dataset that has previously been filtered
    adata = bc.datasets.pbmc3k_filtered()
    ## We get the raw matrix containing all the initial genes, keeping the filtering on the cells
    adata = bc.get_raw(adata)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    /.local/lib/python3.7/site-packages/anndata/compat/__init__.py:161: FutureWarning:

    Moving element from .uns['neighbors']['distances'] to .obsp['distances'].

    This is where adjacency matrices should go now.

    /.local/lib/python3.7/site-packages/anndata/compat/__init__.py:161: FutureWarning:

    Moving element from .uns['neighbors']['connectivities'] to .obsp['connectivities'].

    This is where adjacency matrices should go now.





highly variable gene selection
------------------------------

select highly variable genes (considers correction for gene expression level)


.. code-block:: default


    #define thresholds for highly variable genes
    variable_genes_min_mean = 0.01
    variable_genes_max_mean = 5
    variable_genes_min_disp = 0.4

    #identify genes with variable expression
    filter_result = sc.pp.filter_genes_dispersion(adata.X, min_mean=variable_genes_min_mean, max_mean=variable_genes_max_mean, min_disp=variable_genes_min_disp) 
    sc.pl.filter_genes_dispersion(filter_result)
    nbr_variable_genes = sum(filter_result.gene_subset)
    print('number of variable genes selected ', nbr_variable_genes )

    #perform the actual filtering
    adata = adata[:, filter_result.gene_subset]




.. image:: /auto_examples/preprocessing/images/sphx_glr_plot_pca_neighbors_clustering_001.png
    :alt: plot pca neighbors clustering
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    number of variable genes selected  1897
    /.local/lib/python3.7/site-packages/anndata/_core/anndata.py:1094: FutureWarning:

    is_categorical is deprecated and will be removed in a future version.  Use is_categorical_dtype instead





set random seed
---------------
To get reproducible results you need to define a random seed for all of the stochastic
processes, such as e.g. PCA, neighbors, etc.


.. code-block:: default


    #set random seed
    random_seed = 0








PCA
---


.. code-block:: default


    #log transform our data (is easier to work with numbers like this)
    sc.pp.log1p(adata)

    # Scale data to unit variance and zero mean, and cut-off at max value 10
    sc.pp.scale(adata, max_value=10) 

    #calculate 50 principle components of the dataset
    sc.tl.pca(adata, random_state=random_seed, svd_solver='arpack')

    #visualize the amount of variance explained by each PC
    sc.pl.pca_variance_ratio(adata)

    #visualize the loadings onto the first 3 PCs
    sc.pl.pca_loadings(adata)




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/preprocessing/images/sphx_glr_plot_pca_neighbors_clustering_002.png
          :alt: variance ratio
          :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/preprocessing/images/sphx_glr_plot_pca_neighbors_clustering_003.png
          :alt: PC1, PC2, PC3
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    /.conda/envs/besca_docs/lib/python3.7/site-packages/scanpy/preprocessing/_simple.py:339: UserWarning:

    Revieved a view of an AnnData. Making a copy.





nearest neighbors
-----------------


.. code-block:: default


    sc.pp.neighbors(adata, n_neighbors=15, random_state = random_seed, n_pcs=50)








louvain clustering
------------------


.. code-block:: default


    sc.tl.leiden(adata, random_state=random_seed)








UMAP and t-SNE generation
-------------------------


.. code-block:: default


    #calculate UMAP
    sc.tl.umap(adata, random_state = random_seed)

    #calculate t-SNE
    sc.tl.tsne(adata, random_state = random_seed)








visualize the results
---------------------


.. code-block:: default


    sc.pl.umap(adata, color = ['leiden'])
    sc.pl.tsne(adata, color = ['leiden'])




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/preprocessing/images/sphx_glr_plot_pca_neighbors_clustering_004.png
          :alt: leiden
          :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/preprocessing/images/sphx_glr_plot_pca_neighbors_clustering_005.png
          :alt: leiden
          :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    /.local/lib/python3.7/site-packages/anndata/_core/anndata.py:1192: FutureWarning:

    is_categorical is deprecated and will be removed in a future version.  Use is_categorical_dtype instead






.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 41 minutes  7.887 seconds)


.. _sphx_glr_download_auto_examples_preprocessing_plot_pca_neighbors_clustering.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_pca_neighbors_clustering.py <plot_pca_neighbors_clustering.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_pca_neighbors_clustering.ipynb <plot_pca_neighbors_clustering.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
