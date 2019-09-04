#!/usr/bin/env python
# coding: utf-8

import argparse
import scipy
import anndata

parser = argparse.ArgumentParser(description='Normalisation with scran')
parser.add_argument('-i', '--input_file', required=True)
parser.add_argument('-c', '--cc_markers', help="file containing cell cycle marker genes (assumes that this is Tirosh et al. marker gene file from the scanpy tutorial https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt)") # TODO: split in s genes and g2m genes
parser.add_argument('-o', '--output_file', required=True)
args = parser.parse_args()

import pandas as pd
import scanpy as sc
sc.settings.verbosity = 0
import scIB.preprocessing as pp

adata = sc.read(args.input_file, cache=True)
if isinstance(adata.X, scipy.sparse.csr_matrix):
    adata = anndata.AnnData(X=adata.X.sorted_indices(), obs=adata.obs, var=adata.var)

# Normalisation
print("Normalisation")
pp.normalize(adata)

# Cell cycle scores
if args.cc_markers is not None:
    print("Cell cycle scores")
    s_genes, g2m_genes = pp.cc_tirosh(args.cc_markers)
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

adata.write(args.output_file)
