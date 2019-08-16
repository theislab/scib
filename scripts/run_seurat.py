#!/usr/bin/env python
# coding: utf-8

import argparse

parser = argparse.ArgumentParser(description='Run scanorama for ')
parser.add_argument('-i', '--input_file', required=True)
parser.add_argument('-o', '--output_file', required=True)
parser.add_argument('-b', '--batch', required=True, help='Batch variable')
parser.add_argument('-v', '--hvgs', help='Preselect for HVGs', action='store_true')

args = parser.parse_args()
file = args.input_file
out = args.output_file
batch = args.batch
hvg = args.hvgs


import scanpy as sc
import scIB
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import scipy


# In[3]:


import warnings
warnings.filterwarnings('ignore')
import copy




#batch = 'method'
#hvg = True


#file="/storage/groups/ml01/workspace/group.daniela/Benchmarking_data_integration/merged_adata.h5ad"


adata = sc.read(file)


methods = {}
adatas = {}


# ## Select Highly variable genes
# Set `hvg` to `True`if you want to preselect for less than 4000 HVGs before running the integration methods. This is probably preferable, as many methods work better with a reduced dataset.

if hvg:
    hvgs = scIB.preprocessing.hvg_intersect(adata, batch)
    adata = adata[:,hvgs]


# ## Run the integration method
# The functions for the integration methods are in `scIB.integration`. Generally, the methods expect an anndata object and the batch key as an input. The runtime and memory usage of the functions are meaured using `scIB.metrics.measureTM`. This function returns memory usage in MB, runtime in s and the output of the tested function.


outp = scIB.metrics.measureTM(scIB.integration.runSeurat, adata, batch)


integrated = outp[2][0]


integrated.uns['mem']=outp[0]
integrated.uns['runtime']=outp[1]


# ## Quantifying quality of Integration

# ### Silhouette score


#sil = scIB.metrics.silhouette_score(integrated, batch, 'cell_ontology_class')


# In[15]:


#sc.tl.louvain(adatas['scanorama'], key_added='louvain_post')
#nmi= scIB.metrics.nmi(adatas['scanorama'], 'cell_ontology_class', 'louvain_post')



#pcr_hvg = scIB.metrics.pcr_hvg(adata, adatas['scanorama'], 10, batch)

#pcr_hvg

sc.write(out, integrated)