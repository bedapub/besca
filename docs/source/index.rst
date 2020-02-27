.. _index:

.. toctree::
   :maxdepth: 2
   :caption: contents:
   :hidden:

   besca <besca>
   code examples <auto_examples/index>
   tutorials <tutorials>
   standard pipeline <besca_standard_pipeline>
   adding new functions to besca <adding_new_functions>

Welcome to besca's documentation!
=================================

The besca (BEDA's single cell sequencing analysis) package contains many usefull python functions to use for your single-cell analysis. 

.. image:: _images/besca_outline.jpg
   :align: center

The package has been grouped into 5 categories:

- :ref:`preprocessing functions<preprocessing-functions>`: this submodule contains all functions relevant to data preprocessing  
- :ref:`plotting functions <plotting-functions>`: additional plot types not available in the standard scanpy package  
- :ref:`tools <tools-functions>`: contains additional tools to e.g. perform differential gene analysis
- :ref:`import<import-functions>`/:ref:`export <export-functions>`: collection of functions to export/load data from the FAIR data format
- :ref:`standardworkflow <standardworkflow-functions>`: contains functions optimized for use in our standard scsequencing analysis pipeline

In addition you will find example code and output (including some short tutorials) :doc:`here <auto_examples/index>`, aswell as extensive documentation on :ref:`adding functions to besca <adding-new-functions>` and :ref:`maintaining the package <besca-maintenance>`. 
