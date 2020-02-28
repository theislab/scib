# Benchmarking atlas-level data integration in single-cell genomics

This repository contains code and analysis for the benchmarking study for data integration tools.
The scib python module is in the folder scIB. It can be installed using `pip install -e .` run in the root directory.
R helper functions for R integration methods can be found in the `R` directory.

The `scripts` folder contains scripts for preparing the data, running the methods, postprocessing and calculation of the metrics.

## Running the integration methods
This package allows to run a multitude of single cell data integration methods in both `R` and `python`.
Before running the methods, the data needs to be preprocessed using the `scripts/runPP.py` script. This script calculates highly variable
genes per batch and converts the data to an R Seurat object if needed.
``` python scripts/runPP.py -i INPUT_FILE -o OUTPUT_FILE -b BATCH_VARIABLE [-v #_OF_HVGs]```
Use the `-r` flag to output an R object and `-s` to scale the data.

To run the integration methods the `scripts/runIntegration.py` script is used for methods implemented in python.
``` python scripts/runIntegration.py -i INPUT_FILE -o OUTPUT_FILE -b BATCH_VARIBABLE [-v #_OF_HVGs] ```

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
