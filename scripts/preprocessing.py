#!/usr/bin/env python
# coding: utf-8

import argparse

parser = argparse.ArgumentParser(description='Normalisation with scran')
parser.add_argument('-i', '--input_file', required=True)
parser.add_argument('-c', '--cc_markers', required=True) # TODO: split in s genes and g2m genes
parser.add_argument('-o', '--output_file', required=True)
parser.add_argument('--flavor', default='cell_ranger', help='For HVG computation. Either "cell_ranger" (for log-transformed data) or "seurat" (for no-transformed data).')
parser.add_argument('--n_top_genes', default=4000, help='Number of HVG to be computed.')
args = parser.parse_args()

import pandas as pd
import scanpy as sc
sc.settings.verbosity = 0
import scIB.preprocessing as pp

adata = sc.read(args.input_file, cache=True)

# Normalisation
print("Normalisation")
pp.normalize(adata)

# Visualisation
print("Visualisation")
pp.reduce_data(adata, hvg=True, flavor=args.flavor, n_top_genes=args.n_top_genes,
               pca=True, neighbors=True, umap=True)

# Cell cycle scores
print("Cell cycle scores")
s_genes, g2m_genes = pp.cc_tirosh(cc_markers)
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

adata.write(args.output_file)
