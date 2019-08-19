#!/usr/bin/env python
# coding: utf-8

import argparse

parser = argparse.ArgumentParser(description='Normalisation with scran')
parser.add_argument('-i', '--input_file', required=True)
parser.add_argument('-c', '--cc_markers') # TODO: split in s genes and g2m genes
parser.add_argument('-o', '--output_file', required=True)
args = parser.parse_args()

import pandas as pd
import scanpy as sc
sc.settings.verbosity = 0
import scIB.preprocessing as pp

adata = sc.read(args.input_file, cache=True)

# Normalisation
print("Normalisation")
pp.normalize(adata)

# Cell cycle scores
if args.cc_markers is not None:
    print("Cell cycle scores")
    s_genes, g2m_genes = pp.cc_tirosh(cc_markers)
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

adata.write(args.output_file)
