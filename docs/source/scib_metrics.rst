Metrics
=======

.. currentmodule:: scib.metrics

This package contains all the metrics used for benchmarking scRNA-seq data integration performance.
The metrics can be classified into biological conservation and batch removal metrics.
For a detailed description of the metrics implemented in this package, please see our `publication`_.

.. _publication: https://doi.org/10.1038/s41592-021-01336-8


Biological Conservation Metrics
-------------------------------

Biological conservation metrics quantify either the integrity of cluster-based metrics based on clustering results of
the integration output, or the difference in the feature spaces of integrated and unintegrated data.
Each metric is scaled to a value ranging from 0 to 1 by default, where larger scores represent better conservation of
the biological aspect that the metric addresses.

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

Batch correction metrics values are scaled by default between 0 and 1, in which larger scores represent better batch
removal.

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
The main function is :func:`~scib.metrics.metrics`, that provides all the parameters for the different metrics.

.. code-block:: python

    scib.metrics.metrics(adata, adata_int, ari=True, nmi=True)


The remaining functions call the :func:`~scib.metrics.metrics` for

Furthermore, :func:`~scib.metrics.metrics()` is wrapped by convenience functions with preconfigured subsets of metrics
based on expected computation time:

+ :func:`~scib.metrics.metrics_fast()` only computes metrics that require little preprocessing
+ :func:`~scib.metrics.metrics_slim()` includes all functions of :func:`~scib.metrics.metrics_fast()` and adds clustering-based metrics
+ :func:`~scib.metrics.metrics_all()` includes all metrics


.. autosummary::
    :toctree: api/

    metrics
    metrics_fast
    metrics_slim
    metrics_all


.. raw:: html

    <h2>Auxiliary Functions</h2>

.. autosummary::
    :toctree: api/

    lisi_graph
    pcr
    pc_regression
