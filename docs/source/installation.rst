Installation
============

We recommend working with environments such as Conda or virtualenv, so that python and R dependencies are in one place.
Please also check out `scib-pipeline <https://github.com/theislab/scib-pipeline.git>`_ for ready-to-use environments and
an end-to-end workflow.

Requirements
------------

+ Linux or UNIX system
+ Python >= 3.7
+ R >= 3.6


Installation with pip
---------------------

The ``scib`` python package is available on `PyPI <https://pypi.org/>`_ and can be installed through

.. code-block:: bash

    pip install scib


Alternatively, you can also install the package directly from GitHub directly via

.. code-block:: bash

    pip install git+https://github.com/theislab/scib.git


Additionally, in order to run the R package ``kBET``, you need to install it through R.

.. code-block:: R

    install.packages('remotes')
    remotes::install_github('theislab/kBET')


.. note::

    By default dependencies for integration methods are not installed due to dependency clashes.
    In order to use integration methods, see


Installing additional packages
------------------------------

This package contains code for running integration methods as well as for evaluating their output. However, due to
dependency clashes, ``scib`` is only installed with the packages needed for the metrics. In order to use the integration
wrapper functions, we recommend to work with different environments for different methods, each with their own
installation of ``scib``. You can install optional Python dependencies via pip as follows:

.. code-block::

    pip install scib[bbknn]  # using BBKNN
    pip install scib[scanorama]  # using Scanorama
    pip install scib[bbknn,scanorama]  # Multiple methods in one go

.. note::

    Zsh often doesn't like square brackets. If you are a zsh user, use quotation marks around any statements containing
    square brackets. For example:

    .. code-block:: bash

        pip install 'scib[bbknn]'


The ``setup.cfg`` of the source code for a full list of Python dependencies.
