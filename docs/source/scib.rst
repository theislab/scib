Package: scib
=============

Import scib as:

::

    import scib


Preprocessing
-------------

.. automodapi:: scib.preprocessing

    :no-heading:

    :skip: plot_count_filter
    :skip: plot_scatter
    :skip: readConos
    :skip: readSeurat
    :skip: saveSeurat


Integration
------------

.. automodapi:: scib.integration

    :no-heading:

    :skip: runBBKNN
    :skip: runCombat
    :skip: runMNN
    :skip: runDESC
    :skip: runSaucie
    :skip: runScanorama
    :skip: runScanvi
    :skip: runScGen
    :skip: runScvi
    :skip: runTrVae
    :skip: runTrVaep
    :skip: issparse


Metrics
-------

.. currentmodule:: scib.metrics

.. automodule:: scib.metrics


Biological Conservation Metrics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Biological conservation metrics return a value ranging from 0 to 1, in which larger scores represent better conservation
of the biological aspect that the metric addresses.

.. autosummary::
    :toctree: api/

    hvg_overlap
    silhouette
    isolated_labels
    nmi
    ari
    cell_cycle
    trajectory_conservation
    clisi_graph


Batch Correction Metrics
^^^^^^^^^^^^^^^^^^^^^^^^
Batch correction metrics return a value ranging from 0 to 1, in which larger scores represent better batch removal.

.. autosummary::
    :toctree: api/

    graph_connectivity
    silhouette_batch
    pcr_comparison
    kBET
    ilisi_graph


Metrics Wrapper Functions
^^^^^^^^^^^^^^^^^^^^^^^^^

For convenience, ``scib`` provides wrapper functions that, given integrated and unintegrated adata objects, apply
multiple metrics and return all the results in a ``pandas.Dataframe``.
The main function is :func:`scib.metrics.metrics`, that provides all the parameters for the different metrics.
The remaining functions call the :func:`scib.metrics.metrics` for preconfigured subsets of metrics based on runtime.

.. autosummary::
    :toctree: api/

    metrics
    metrics_fast
    metrics_slim
    metrics_all