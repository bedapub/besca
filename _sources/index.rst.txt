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

The package has grouped into 6 modules:

- :ref:`preprocessing functions<preprocessing-functions>`: the `pp` module contains functions relevant for data preprocessing.
- :ref:`plotting functions <plotting-functions>`: the `pl` module offers additional plot types not available in the scanpy package  
- :ref:`tools <tools-functions>`: the `tl` module contains additional tools, for instance tools to perform differential gene analysis
- :ref:`standardworkflow <standardworkflow-functions>`: the `st` module contains functions optimized for besca standard single-cell sequencing analysis pipeline
- :ref:`import<import-functions>`/:ref:`export <export-functions>`: the `Import` and `export` modules are collection of functions to export/load data from the FAIR data format

In addition you will find example code and output (including some short tutorials) :doc:`here <auto_examples/index>`, as well as extensive documentation on :ref:`adding functions to besca <adding-new-functions>`.
