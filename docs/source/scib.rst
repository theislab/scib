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


Biological Conservation Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


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


Batch Correction Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^


.. autosummary::
    :toctree: api/

    graph_connectivity
    silhouette_batch
    pcr_comparison
    kBET
    ilisi_graph


Metrics Wrapper Functions
^^^^^^^^^^^^^^^^^^^^^^^^^

Wrapper functions are available that compute multiple metrics and return them as a ``pandas.Dataframe``.

.. autosummary::
    :toctree: api/

    metrics
    metrics_fast
    metrics_slim
    metrics_all