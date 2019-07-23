#!/bin/env python

### D. C. Strobl, M. Müller; 2019-07-23

""" This module provides a toolkit for running a large range of single cell data integration methods
    as well as tools and metrics to benchmark them.
"""

import scanpy as sc
import anndata
import numpy as np

# first, some checker functions for data sanity
def checkAdata(adata):
    if type(adata) is not anndata.AnnData:
        raise TypeError('Input is not a valid AnnData object')

def checkBatch(batch, obs):
    if batch not in obs:
        raise ValueError('Selected batch column is not in obs')
    else:
        nBatch = obs[batch].nunique()
        print('Object contains '+str(nBatch)+' batches.')

def splitBatches(adata, batch):
    split = []
    for i in adata.obs[batch].unique():
        split.append(adata[adata.obs[batch]==i])
    return split


# functions for running the methods

def runScanorama(adata, batch, hvgPre):
    import scanorama
    checkAdata(adata)
    checkBatch(batch, adata.obs)
    split = splitBatches(adata.copy(), batch)
    emb, corrected = scanorama.correct_scanpy(split, return_dimred=True)
    corrected = corrected[0].concatenate(corrected[1:])
    emb = np.concatenate(emb, axis=0)

    return emb, corrected


if __name__=="__main__":
    adata = sc.read('testing.h5ad')
    emb, corrected = runScanorama(adata, 'method', False)
    print(emb)
    print(corrected)


        

