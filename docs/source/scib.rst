Package: scib
=============

The package contains several modules for preprocessing an ``anndata`` object, running integration methods and
evaluating the resulting.
For preprocessing, ``scib.preprocessing`` (or ``scib.pp``) contains functions for normalising, scaling or batch-aware
selection of highly variable genes.
Functions for the integration methods are in ``scib.integration`` or for short ``scib.ig`` and metrics are under
``scib.metrics`` (or ``scib.me``).


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

Integration method functions require the preprocessed ``anndata`` object (here ``adata``) and the name of the batch column
in ``adata.obs`` (here ``'batch'``).
The methods can be called using the following, where ``<method>`` is the name of the integration method.

.. code-block:: python

    scib.ig.<method>(adata, batch='batch')


For example, in order to run Scanorama, on a dataset, call:

.. code-block:: python

    scib.ig.scanorama(adata, batch='batch')

.. warning::

    The following notation is deprecated.

    .. code-block:: python

        scib.integration.run<method>(adata, batch='batch')

    Please use the snake_case naming without the ``run`` prefix.

Some integration methods (``scgen``, ``scanvi``) also use cell type labels as input. For these, you need to additionally
provide the corresponding label column of ``adata.obs`` (here ``cell_type``).

.. code-block:: python

    scib.ig.scgen(adata, batch='batch', cell_type ='cell_type')
    scib.ig.scanvi(adata, batch='batch', labels ='cell_type')


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