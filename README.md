# Benchmarking atlas-level data integration in single-cell genomics

This repository contains the code for the `scib` package used in our benchmarking study for data integration tools.
In [our study](https://doi.org/10.1038/s41592-021-01336-8), we benchmark 16 methods (see Tools) with 4 combinations of
preprocessing steps leading to 68 methods combinations on 85 batches of gene expression and chromatin accessibility data.

![Workflow](https://raw.githubusercontent.com/theislab/scib/main/figure.png)

## Resources

+ The git repository of the [`scib` package](https://github.com/theislab/scib) and its [documentation](https://scib.readthedocs.io/).
+ The reusable pipeline we used in the study can be found in the
  separate [scib pipeline](https://github.com/theislab/scib-pipeline.git) repository. It is reproducible and automates
  the computation of preprocesssing combinations, integration methods and benchmarking metrics.
+ On our [website](https://theislab.github.io/scib-reproducibility) we visualise the results of the study.
+ For reproducibility and visualisation we have a dedicated
  repository: [scib-reproducibility](https://github.com/theislab/scib-reproducibility).

### Please cite:

Luecken, M.D., Büttner, M., Chaichoompu, K. et al. Benchmarking atlas-level data integration in single-cell genomics.
Nat Methods 19, 41–50 (2022). [https://doi.org/10.1038/s41592-021-01336-8](https://doi.org/10.1038/s41592-021-01336-8)

## `scib` Package

We created the python package called `scib` that uses `scanpy` to streamline the integration of single-cell datasets
and evaluate the results.
For evaluating the integration quality it provides a number of metrics.

The `scib` python package is available on [PyPI](https://pypi.org/) and can be installed through

```
pip install scib
```

## Metrics

For a detailed description of the metrics implemented in this package, please see
our [publication](https://doi.org/10.1038/s41592-021-01336-8).

### Batch removal metrics include:

- Principal component regression `scib.metrics.pcr_comparison()`
- Batch ASW `scib.metrics.silhouette_batch()`
- K-nearest neighbour batch effect `scib.metrics.kBET()`
- Graph connectivity `scib.metrics.graph_connectivity()`
- Graph iLISI `scib.metrics.ilisi_graph()`

### Biological conservation metrics include:

- Normalised mutual information `scib.metrics.nmi()`
- Adjusted Rand Index `scib.metrics.ari()`
- Cell type ASW `scib.metrics.silhouette()`
- Isolated label score F1 `scib.metrics.isolated_labels()`
- Isolated label score ASW `scib.metrics.isolated_labels()`
- Cell cycle conservation `scib.metrics.cell_cycle()`
- Highly variable gene conservation `scib.metrics.hvg_overlap()`
- Trajectory conservation `scib.metrics.trajectory_conservation()`
- Graph cLISI `scib.metrics.clisi_graph()`

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