[![Stars](https://img.shields.io/github/stars/theislab/scib?logo=GitHub&color=yellow)](https://github.com/theislab/scib/stargazers)
[![PyPI](https://img.shields.io/pypi/v/scib?logo=PyPI)](https://pypi.org/project/scib)
[![PyPIDownloads](https://pepy.tech/badge/scib)](https://pepy.tech/project/scib)
[![Build Status](https://github.com/theislab/scib/actions/workflows/test.yml/badge.svg)](https://github.com/theislab/scib/actions/workflows/test.yml)
[![Documentation](https://readthedocs.org/projects/scib/badge/?version=latest)](https://scib.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/theislab/scib/branch/main/graph/badge.svg?token=M1nuTpAxyS)](https://codecov.io/gh/theislab/scib)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)

# Benchmarking atlas-level data integration in single-cell genomics

This repository contains the code for the `scib` package used in our benchmarking study for data integration tools.
In [our study](https://doi.org/10.1038/s41592-021-01336-8), we benchmark 16 methods (see Tools) with 4 combinations of
preprocessing steps leading to 68 methods combinations on 85 batches of gene expression and chromatin accessibility
data.

![Workflow](https://raw.githubusercontent.com/theislab/scib/main/figure.png)

## Resources

- The git repository of the [`scib` package](https://github.com/theislab/scib) and
  its [documentation](https://scib.readthedocs.io/).
- The reusable pipeline we used in the study can be found in the
  separate [scib pipeline](https://github.com/theislab/scib-pipeline.git) repository. It is reproducible and automates
  the computation of preprocesssing combinations, integration methods and benchmarking metrics.
- On our [website](https://theislab.github.io/scib-reproducibility) we visualise the results of the study.
- For reproducibility and visualisation we have a dedicated
  repository: [scib-reproducibility](https://github.com/theislab/scib-reproducibility).

### Please cite:

Luecken, M.D., Büttner, M., Chaichoompu, K. et al. Benchmarking atlas-level data integration in single-cell genomics.
Nat Methods 19, 41–50 (2022). [https://doi.org/10.1038/s41592-021-01336-8](https://doi.org/10.1038/s41592-021-01336-8)

## Package: scib

We created the python package called `scib` that uses `scanpy` to streamline the integration of single-cell datasets and
evaluate the results. The package contains several modules for preprocessing an `anndata` object, running integration
methods and evaluating the resulting using a number of metrics. For preprocessing, `scib.preprocessing` (or `scib.pp`)
contains functions for normalising, scaling or batch-aware selection of highly variable genes. Functions for the
integration methods are in `scib.integration` or for short `scib.ig` and metrics are under
`scib.metrics` (or `scib.me`).

The `scib` python package is available on [PyPI](https://pypi.org/) and can be installed through

```commandline
pip install scib
```

or from the source code.
```commandline
pip install .
```

Import `scib` in python:

```python
import scib
```

### Optional Dependencies

The package contains optional dependencies that need to be installed manually if needed.
These include R dependencies (`rpy2`, `anndata2ri`) which require an installation of R integration method packages.
All optional dependencies are listed under `setup.cfg` under `[options.extras_require]` and can be installed through pip.

e.g. for installing `rpy2` and `bbknn` dependencies:
```commandline
pip install 'scib[rpy2,bbknn]'
```

Optional dependencies outside of python need to be installed separately.
For instance, in order to run kBET, install it via the following command in R:

```R
remotes::install_github('theislab/kBET')
```

## Metrics

We implemented different metrics for evaluating batch correction and biological conservation in the `scib.metrics`
module.

<table class="docutils align-default">
  <colgroup>
    <col style="width: 50%" />
    <col style="width: 50%" />
  </colgroup>
  <thead>
    <tr class="row-odd"><th class="head"><p>Biological Conservation</p></th>
      <th class="head"><p>Batch Correction</p></th>
    </tr>
  </thead>
  <tbody>
    <tr class="row-even" >
      <td><ul class="simple">
        <li><p>Cell type ASW</p></li>
        <li><p>Cell cycle conservation</p></li>
        <li><p>Graph cLISI</p></li>
        <li><p>Adjusted rand index (ARI) for cell label</p></li>
        <li><p>Normalised mutual information (NMI) for cell label</p></li>
        <li><p>Highly variable gene conservation</p></li>
        <li><p>Isolated label ASW</p></li>
        <li><p>Isolated label F1</p></li>
        <li><p>Trajectory conservation</p></li>
      </ul></td>
      <td><ul class="simple">
        <li><p>Batch ASW</p></li>
        <li><p>Principal component regression</p></li>
        <li><p>Graph iLISI</p></li>
        <li><p>Graph connectivity</p></li>
        <li><p>kBET (K-nearest neighbour batch effect)</p></li>
      </ul></td>
    </tr>
  </tbody>
</table>

For a detailed description of the metrics implemented in this package, please see our
[publication](https://doi.org/10.1038/s41592-021-01336-8) and the package [documentation](https://scib.readthedocs.io/).

## Integration Tools

Tools that are compared include:

- [BBKNN](https://github.com/Teichlab/bbknn) 1.3.9
- [Combat](https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.combat.html) [paper](https://academic.oup.com/biostatistics/article/8/1/118/252073)
- [Conos](https://github.com/hms-dbmi/conos) 1.3.0
- [DESC](https://github.com/eleozzr/desc) 2.0.3
- [FastMNN](https://bioconductor.org/packages/batchelor/) (batchelor 1.4.0)
- [Harmony](https://github.com/immunogenomics/harmony) 1.0
- [LIGER](https://github.com/MacoskoLab/liger) 0.5.0
- [MNN](https://github.com/chriscainx/mnnpy) 0.1.9.5
- [SAUCIE](https://github.com/KrishnaswamyLab/SAUCIE)
- [Scanorama](https://github.com/brianhie/scanorama) 1.7.0
- [scANVI](https://github.com/chenlingantelope/HarmonizationSCANVI) (scVI 0.6.7)
- [scGen](https://github.com/theislab/scgen) 1.1.5
- [scVI](https://github.com/YosefLab/scVI) 0.6.7
- [Seurat v3](https://github.com/satijalab/seurat) 3.2.0 CCA (default) and RPCA
- [TrVae](https://github.com/theislab/trvae) 0.0.1
- [TrVaep](https://github.com/theislab/trvaep) 0.1.0

## Development

For developing this package, please make sure to install additional dependencies so that you can use `pytest` and
`pre-commit`.

```shell
pip install -e '.[test,dev]'
```

Please refer to the `setup.cfg` for more optional dependencies.

Install `pre-commit` to the repository for running it automatically every time you commit in git.

```shell
pre-commit install
```
