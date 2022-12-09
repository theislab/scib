API
===

Preprocessing
-------------

.. currentmodule:: scib.preprocessing

Preprocessing functions are relevant both for preparing the data for integration as well as postprocessing the
integration output.

The most relevant preprocessing steps are:

+ Normalization
+ Scaling, batch-aware
+ Highly variable gene selection, batch-aware
+ Cell cycle scoring
+ Principle component analysis (PCA)
+ k-nearest neighbor graph (kNN graph)
+ UMAP
+ Clustering

Note that some preprocessing steps depend on each other.
Please refer to the `best_practices`_ for more details.

.. _best_practices: sc-best-practices.org

.. autosummary::
    :toctree: api/

    normalize
    scale_batch
    hvg_intersect
    hvg_batch
    score_cell_cycle
    reduce_data


Integration
-----------

Integration method functions require the preprocessed ``anndata`` object (here ``adata``) and the name of the batch column
in ``adata.obs`` (here ``'batch'``).
The methods can be called using the following, where ``integration_method`` is the name of the integration method.

.. code-block:: python

    scib.ig.integration_method(adata, batch="batch")


For example, in order to run Scanorama, on a dataset, call:

.. code-block:: python

    scib.ig.scanorama(adata, batch="batch")

.. warning::

    The following notation is deprecated.

    .. code-block:: python

        scib.integration.runIntegrationMethod(adata, batch="batch")

    Please use the snake_case naming without the ``run`` prefix.

Some integration methods (e.g. :func:`~scib.integration.scgen`, :func:`~scib.integration.scanvi`) also use cell type
labels as input.
For these, you need to additionally provide the corresponding label column of ``adata.obs`` (here ``cell_type``).

.. code-block:: python

    scib.ig.scgen(adata, batch="batch", cell_type="cell_type")
    scib.ig.scanvi(adata, batch="batch", labels="cell_type")


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

This package contains all the metrics used for benchmarking scRNA-seq data integration performance.
They can be applied on the integrated as well as the unintegrated data and  can be classified into biological
conservation and batch removal metrics.
For a detailed description of the metrics implemented in this package, please see our `publication`_.

.. _publication: https://doi.org/10.1038/s41592-021-01336-8

Most metrics require specific inputs that need to be preprocessed, which is described in detail under :ref:`preprocessing`.


Biological Conservation Metrics
```````````````````````````````

Biological conservation metrics quantify either the integrity of cluster-based metrics based on clustering results of
the integration output, or the difference in the feature spaces of integrated and unintegrated data.
Each metric is scaled to a value ranging from 0 to 1 by default, where larger scores represent better conservation of
the biological aspect that the metric addresses.

.. autosummary::
    :toctree: api/

    ari
    cell_cycle
    clisi_graph
    hvg_overlap
    isolated_labels_asw
    isolated_labels_f1
    nmi
    silhouette
    trajectory_conservation


Batch Correction Metrics
````````````````````````

Batch correction metrics values are scaled by default between 0 and 1, in which larger scores represent better batch
removal.

.. autosummary::
    :toctree: api/

    graph_connectivity
    ilisi_graph
    kBET
    pcr_comparison
    silhouette_batch


Metrics Wrapper Functions
`````````````````````````

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


Auxiliary Functions
```````````````````

Some parts of metrics can be used individually, these are listed below.

.. autosummary::
    :toctree: api/

    lisi_graph
    pcr
    pc_regression
