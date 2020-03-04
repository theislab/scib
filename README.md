# Benchmarking atlas-level data integration in single-cell genomics

This repository contains code and analysis for the benchmarking study for data integration tools.
The scib python module is in the folder scIB. It can be installed using `pip install -e .` run in the root directory.
R helper functions for R integration methods can be found in the `R` directory.

## Installation
To reproduce the results from this study, three different conda environments are needed.
This is due to package incompatibilities from different integration methods and other necessary packages.
There are different environments for the python integration methods, the R integration methods and
the conversion of R data types to anndata objects.

For the installation of conda, follow [these](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) instructions
or use your system's package manager. The environments have only been tested on linux operating systems
although it might be possible to run the pipeline using Mac OS.

To create the conda environments use the `.yml` files in the `envs` directory.
To install the envs, use `conda env create -f FILENAME.yml`.

The `scripts` folder contains scripts for preparing the data, running the methods, postprocessing and calculation of the metrics.

## Running the integration methods
This package allows to run a multitude of single cell data integration methods in both `R` and `python`.
We use [Snakemake](https://snakemake.readthedocs.io/en/stable/) to run the pipeline.
The parameters of the run are configured using the `config.yaml` file.
See the `DATA_SCENARIOS` section to change the data used for integration.
The script expects one `.h5ad` file containing all batches per data scenario.

To load the config file run `snakemake --configfile config.yaml`.
Define the number of CPU threads you want to use with `snakemake --cores N_CORES`. To produce an overview of tasks that will be run, use `snakemake -n`.
To run the pipeline, simply run `snakemake`.
## Tools
Tools to be compared include:
- Seurat v2
- [Seurat v3](https://github.com/satijalab/seurat)
- [TrVae](https://github.com/theislab/trvae)
- [scVI](https://github.com/YosefLab/scVI)
- [scANVI](https://github.com/chenlingantelope/HarmonizationSCANVI)
- [CONOS](https://github.com/hms-dbmi/conos) [tutorial](https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/conos.html)
- [MNN](https://github.com/chriscainx/mnnpy)
- [Scanorama](https://github.com/brianhie/scanorama)
- RISC
- [LIGER](https://github.com/MacoskoLab/liger)
- [BBKNN](https://github.com/Teichlab/bbknn)
- [Harmony](https://github.com/immunogenomics/harmony)
- [scMerge](https://github.com/SydneyBioX/scMerge)
- [scAlign](https://github.com/quon-titative-biology/scAlign)
- BBKNN + [scAEspy](https://gitlab.com/cvejic-group/scaespy)?


## Data
Data scenarios to be compared are:

1. Pancreas data

2. Mouse brain data

3. Immune cell data across species and perturbations

4. Mouse atlases:
- Tabula Muris
- Mouse atlas

5. Mouse brain scATAC-seq data

Possible extensions:
- potentially 1 more scATAC-seq data scenario
- If possible, we will attempt to merge scATAC and scRNA-seq datasets in mouse brain as well.
- lung atlas mapping using our published datasets (inter-individual variation)

## Metrics
+ silhouette score
+ kBET
+ PCA regression
+ normalized mutual information NMI [github](https://github.com/aaronmcdaid/Overlapping-NMI)

## Strategy

1. Datasets are normalized and pre-processed separately. QC from the publications will be used wherever possible.
2. Cell labels must be harmonized between datasets (ideally on the basis of the cell ontology). This can be done by matching marker gene sets (overlap metrics).
3. Data integration is performed with default settings
4. Integrated datasets are scored based on metrics.
