Integration
===========

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
