Metrics
=======

.. currentmodule:: scib.metrics

This package contains all the metrics used for benchmarking scRNA-seq data integration performance.
The metrics can be classified into biological conservation and batch removal metrics.
For a detailed description of the metrics implemented in this package, please see our `publication`_.

.. _publication: https://doi.org/10.1038/s41592-021-01336-8

Preprocessing data for metrics
------------------------------

The metrics can either be used to evaluate a dataset or to evaluate the integration performance on different
representations of the data.
Data representations can be one of the following:

1. **Feature space:** a feature-level expression matrix (usually normalised & log-transformed counts or corrected counts from integration)
2. **Embedding space:** PCA of the feature-level expression matrix or an embedding from an integration method
3. **kNN graph space:** a kNN graph on a PCA or an embedding

Unintegrated feature matrices can simply be reduced to a PCA embedding or kNN graph using the according ``scanpy``
functions or the  wrapper :func:`~scib.preprocessing.reduce_data`.
The integrated output can come in one of the three representations, which are denoted as "native representation" in the
figure below.
Each metric assumes that the data representation it requires is available in the ``anndata`` object.
If that is not the case, the data needs to be preprocessed accordingly or, for cases where that is not possible, the
metric cannot be applied to that data representation.

This means that native representations in kNN graph space can only be evaluated by graph based metrics, while embedding
space outputs can be evaluated by embedding and graph based metrics, the latter requiring processing for the kNN graph.
From a feature space native representation one can compute embedding space (PCA on expression matrix, optionally with
highly variable gene selection) and kNN graph (on PCA) and all metrics are applicable.
These relationships are visualised in the figure below.

.. figure:: _static/metrics_workflow.png
   :alt: metrics workflow

   **Overview of metrics and processing steps for different data representations.**
   Metrics are arranged depending on the data representation they require.
   The different data representations are separated by dotted lines in the figure.
   Native data representations (from integration output) are depicted in green boxes, while the processed
   representations are in white boxes.
   The preprocessing commands are denoted on the arrows between data representations, highly variable gene selection is
   not shown here explicitly.
   Metrics that are evaluate differently for different data representations are depicted by metrics on top of the dotted
   lines.

The data representation not only determines which metrics can be evaluated, but also how they are computed.
For instance, the principle component regression :func:`~scib.metrics.pcr_comparison` requires a PCA and variance
contributions, which would classically be computed a full feature matrix but can also be computed on an integrated
embedding matrix.


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
    isolated_labels_f1
    isolated_labels_asw
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
