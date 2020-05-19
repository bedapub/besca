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

The besca (BEDA's single-cell sequencing analysis) package builds upon the scanpy library and offers additional data structures and functions for single-cell analysis.

.. image:: _images/besca_outline.jpg
   :align: center

The package has grouped into 5 sub-modules:

- :ref:`preprocessing functions<preprocessing-functions>`: this sub-module contains all functions relevant to data preprocessing  
- :ref:`plotting functions <plotting-functions>`: additional plot types not available in the standard scanpy package  
- :ref:`tools <tools-functions>`: contains additional tools to e.g. perform differential gene analysis
- :ref:`import<import-functions>`/:ref:`export <export-functions>`: collection of functions to export/load data from the FAIR data format
- :ref:`standardworkflow <standardworkflow-functions>`: contains functions optimized for use in our standard single-cell sequencing analysis pipeline

In addition you will find example code and output (including some short tutorials) :doc:`here <auto_examples/index>`, as well as extensive documentation on :ref:`adding functions to besca <adding-new-functions>` and :ref:`maintaining the package <besca-maintenance>`.
