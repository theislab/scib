#!/usr/bin/env python
# coding: utf-8

import argparse
import scipy
import anndata

parser = argparse.ArgumentParser(description='Normalisation with scran')
parser.add_argument('-i', '--input_file', required=True)
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

adata.write(args.output_file)
