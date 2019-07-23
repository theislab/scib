# Benchmarking_data_integration

This repository contains code and analysis for the benchmarking study for data integration tools.

Tools to be compared include:
- Seurat v2
- [Seurat v3](https://github.com/satijalab/seurat)
- [scGen](https://github.com/theislab/scgen)
- [scVI](https://github.com/YosefLab/scVI)
- [CONOS](https://github.com/hms-dbmi/conos)
- [MNN](https://github.com/chriscainx/mnnpy)
- [Scanorama](https://github.com/brianhie/scanorama)
- RISC
- [LIGER](https://github.com/MacoskoLab/liger)
- [BBKNN](https://github.com/Teichlab/bbknn)
- [Harmony](https://github.com/immunogenomics/harmony)
- [scMerge](https://github.com/SydneyBioX/scMerge)


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


Strategy:

1. Datasets are normalized and pre-processed separately. QC from the publications will be used wherever possible.
2. Cell labels must be harmonized between datasets (ideally on the basis of the cell ontology). This can be done by matching marker gene sets (overlap metrics).
3. Data integration is performed with default settings
4. Integrated datasets are scored based on metrics.
