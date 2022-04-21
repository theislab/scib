Metrics
=======

.. currentmodule:: scib.metrics

.. automodule:: scib.metrics


Biological Conservation Metrics
-------------------------------

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
------------------------

Batch correction metrics return a value ranging from 0 to 1, in which larger scores represent better batch removal.

.. autosummary::
    :toctree: api/

    graph_connectivity
    silhouette_batch
    pcr_comparison
    kBET
    ilisi_graph


Metrics Wrapper Functions
-------------------------

For convenience, ``scib`` provides wrapper functions that, given integrated and unintegrated adata objects, apply
multiple metrics and return all the results in a ``pandas.Dataframe``.
The main function is :func:`scib.metrics.metrics`, that provides all the parameters for the different metrics.

.. code-block:: python

    scib.metrics.metrics(adata, adata_int, ari=True, nmi=True)


The remaining functions call the :func:`scib.metrics.metrics` for

Furthermore, :func:`scib.metrics.metrics()` is wrapped by convenience functions with preconfigured subsets of metrics
based on expected computation time:

+ :func:`scib.me.metrics_fast()` only computes metrics that require little preprocessing
+ :func:`scib.me.metrics_slim()` includes all functions of :func:`scib.me.metrics_fast()` and adds clustering-based metrics
+ :func:`scib.me.metrics_all()` includes all metrics


.. autosummary::
    :toctree: api/

    metrics
    metrics_fast
    metrics_slim
    metrics_all