# Benchmarking atlas-level data integration in single-cell genomics

This repository contains the code for our benchmarking study for data integration tools.
In [our study](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v1), we benchmark 16 methods ([see here](##Tools)) with 4 combinations of preprocessing steps leading to 68 
methods combinations on 85 batches of gene expression and chromatin accessibility data.

![Workflow](./figure.png)

## Package: `scIB`
We created the python package called `scIB` that uses `scanpy` to streamline the integration of single-cell datasets and evaluate the results.
The evaluation of integration quality is based on a number of metrics.

### Installation
The `scIB` python package is in the folder scIB.
You can install it from the root of this repository using

```
pip install .
```

Additionally, in order to run the R package `kBET`, you need to install it through R.

```R
devtools::install_github('theislab/kBET')
```

We recommend to use a conda environment or something similar, so that python and R dependencies are in one place.
Please also check out [scIB pipeline](https://github.com/theislab/scib-pipeline.git) for ready-to-use environments.

### Package structure
The package contains several modules for the different steps of the integration and benchmarking pipeline.
Functions for the integration methods are in `scIB.integration`. The methods are called using `scIB.integration.runMETHOD(adata, BATCH)`.
`scIB.preprocessing` contains methods for preprocessing of the data such as normalisation, scaling or highly variable gene selection per batch.
The metrics are located at `scIB.metrics`. To run multiple metrics in one run, use the `scIB.metrics.metrics()` function.

### Metrics
For a detailed description of the metrics implemented in this package, please see the [manuscript](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2).

#### Batch removal metrics include:

- Principal component regression (`pcr_comparison()`)
- Batch ASW (`silhouette()`)
- K-nearest neighbour batch effect (`kBET()`)
- Graph connectivity (`graph_connectivity()`)
- Graph iLISI (`lisi_graph()`)

#### Biological conservation metrics include:

- Normalised mutual information (`nmi()`)
- Adjusted Rand Index (`ari()`)
- Cell type ASW (`silhouette_batch()`)
- Isolated label score F1 (`isolated_labels()`)
- Isolated label score ASW (`isolated_labels()`)
- Cell cycle conservation (`cell_cycle()`)
- Highly variable gene conservation (`hvg_overlap()`)
- Trajectory conservation (`trajectory_conservation()`)
- Graph cLISI (`lisi_graph()`)

## Tools
Tools that are compared include:
- [Scanorama](https://github.com/brianhie/scanorama)
- [scANVI](https://github.com/chenlingantelope/HarmonizationSCANVI)
- [FastMNN](https://bioconductor.org/packages/batchelor/)
- [scGen](https://github.com/theislab/scgen)
- [BBKNN](https://github.com/Teichlab/bbknn)
- [scVI](https://github.com/YosefLab/scVI)
- [Seurat v3 (CCA and RPCA)](https://github.com/satijalab/seurat)
- [Harmony](https://github.com/immunogenomics/harmony)
- [Conos](https://github.com/hms-dbmi/conos)
- [Combat](https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.combat.html) [paper](https://academic.oup.com/biostatistics/article/8/1/118/252073)
- [MNN](https://github.com/chriscainx/mnnpy)
- [TrVae](https://github.com/theislab/trvae)
- [DESC](https://github.com/eleozzr/desc)
- [LIGER](https://github.com/MacoskoLab/liger)
- [SAUCIE](https://github.com/KrishnaswamyLab/SAUCIE)

## Resources

+ On our [website](https://theislab.github.io/scib-reproducibility) we visualise the results of the study.

+ The reusable pipeline we used in the study can be found in the separate [scIB pipeline](https://github.com/theislab/scib-pipeline.git) repository. It is reproducible and automates the computation of preprocesssing combinations, integration methods and benchmarking metrics.

+ For reproducibility and visualisation we have a dedicated repository: [scib-reproducibility](https://github.com/theislab/scib-reproducibility).

### Please cite:  

_**Benchmarking atlas-level data integration in single-cell genomics.**  
MD Luecken, M Büttner, K Chaichoompu, A Danese, M Interlandi, MF Mueller, DC Strobl, L Zappia, M Dugas, M Colomé-Tatché, FJ Theis
bioRxiv 2020.05.22.111161; doi: https://doi.org/10.1101/2020.05.22.111161_ 
