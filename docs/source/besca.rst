besca
=====

.. _helper-functions:

helper functions
----------------

.. autosummary::
	:toctree: helper_functions

	besca.get_raw
	besca.subset_adata
	besca.convert_ensembl_to_symbol
	besca.convert_symbol_to_ensembl
	besca.get_raw
	besca.get_means
	besca.get_gmeans
	besca.concate_adata

.. _preprocessing-functions:

preprocessing
-------------
.. automodsumm:: besca.pp
	:toctree: preprocessing

.. _plotting-functions:

plotting
--------
.. automodsumm:: besca.pl
	:toctree: plotting

.. _tools-functions:

tools
-----
.. automodsumm:: besca.tl
	:toctree: tools

.. _toolkits-functions:

toolkits
^^^^^^^^

batch correction
++++++++++++++++
Collection of functions to perform batch correction.

.. automodsumm:: besca.tl.bcor
	:toctree: bcor

differential gene expression
++++++++++++++++++++++++++++
Collection of functions to aid in differential gene expression analysis.

.. automodsumm:: besca.tl.dge
	:toctree:  dge

signature scoring
+++++++++++++++++
Collection of functions to aid in signature scoring.

.. automodsumm:: besca.tl.sig
        :toctree:  sig

reclustering
++++++++++++
Collection of functions to perform reclustering on selected subclusters.

.. automodsumm:: besca.tl.rc
	:toctree: reclustering

.. _import-functions:

auto-annot
++++++++++++
Collection of functions to perform auto-annot : annotating a sc datasets based on a reference one.

.. automodsumm:: besca.tl.auto_annot
	:toctree: auto_annot


Import
------
.. automodsumm:: besca.Import
	:toctree: import

.. _export-functions:

export
------
.. automodsumm:: besca.export
	:toctree: export

.. _standardworkflow-functions:

standardworkflow
----------------
.. automodsumm:: besca.st
	:toctree: standardworkflow
