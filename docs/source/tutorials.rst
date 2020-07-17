.. _tutorials:

tutorials
=========


.. toctree::
   :hidden:
   :maxdepth: 2
   
   tutorials/notebook1_data_processing.ipynb
   tutorials/notebook2_celltype_annotation.ipynb
   tutorials/notebook3_batch_correction.ipynb
   tutorials/auto_annot_tutorial.ipynb
   tutorials/bescape_tutorial.ipynb
   tutorials/adata_to_eset.ipynb
   

single cell sequencing general tutorials
---------------------------------------------------------------
   
This tutorial will give you a general introduction into single-cell sequencing analysis using scanpy and besca. It is divided into three seperate notebooks that should be looked at in consecutive order (they build up on results from the previous notebooks).

**Part 1**: :doc:`data processing <tutorials/notebook1_data_processing>` 

**Part 2**: :doc:`celltype annotation <tutorials/notebook2_celltype_annotation>`

**Part 3**: :doc:`batch correction <tutorials/notebook3_batch_correction>`

After looking at the introductory tutorial you can download the hands-on tutorial yourself from `here <github.com/bedapub/besca/blob/master/docs/source/tutorials/scRNAseq_tutorial.ipynb>`_  (please save the link as a :code:`.ipynb` file) and compare with the results published :doc:`here <tutorials/scRNAseq_tutorial>`.


single cell auto_annot tutorial for cell type annotation
---------------------------------------------------------------

We also provide a tutorial for the auto_annot package, which allows to automatically annotate cell types useing supervised machine learning, you can download it from `here <github.com/bedapub/besca/blob/master/docs/source/tutorials/auto_annot_tutorial.ipynb>`_  (please save the link as a :code:`.ipynb` file) and compare with the results published :doc:`here <tutorials/auto_annot_tutorial>`.


Bescape: Cell deconvolution tutorial
---------------------------------------------------------------
Bescape (BESCA proportion estimator) is a deconvolution module. It utilises single-cell annotations coming from the BESCA workflow to build a Gene Expression Profile (GEP). This GEP is used as a basis vector to deconvolute bulk RNA samples i.e. predict cell type proportions within a sample.

Please find the deconvolution tutorial :doc:`here <tutorials/bescape_tutorial>`.

Some deconvolution methods provided by Bescape are written in R. Thus, we need to convert the AnnData objects to R ExpressionSet objects. This has been semi-automated in the following tutorial: :doc:`Converting AnnData to ExpressionSet <tutorials/adata_to_eset>`

