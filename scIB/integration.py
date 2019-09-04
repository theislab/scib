#!/bin/env python

### D. C. Strobl, M. Müller; 2019-07-23

""" This module provides a toolkit for running a large range of single cell data integration methods
    as well as tools and metrics to benchmark them.
"""

import scanpy as sc
import scipy as sp
#import numpy as np
from scIB.utils import *
from memory_profiler import profile
import os
import pandas as pd
import anndata

import rpy2.rinterface_lib.callbacks
import logging
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR) # Ignore R warning messages
import rpy2.robjects as ro
import anndata2ri

# functions for running the methods

def runScanorama(adata, batch, hvg = None):
    import scanorama
    checkSanity(adata, batch, hvg)
    split = splitBatches(adata.copy(), batch)
    emb, corrected = scanorama.correct_scanpy(split, return_dimred=True)
    corrected = corrected[0].concatenate(corrected[1:])
    emb = np.concatenate(emb, axis=0)
    corrected.obsm['X_pca']= emb
    corrected.uns['emb']=True

    return corrected

def runScGen(adata, batch, cell_type, n_top_genes=4000, model_path='./models/batch', epochs=100, hvg=None):
    checkSanity(adata, batch, hvg)
    import scgen
    
    ### Edited by Kris, please verify. 04.09.2019 17:45
    
    sc.tl.louvain(adata, resolution=0.5)
    
    ### End edited by Kris
    
    # save cell_types for later
    cell_types = adata.obs[cell_type].copy()
    batches = adata.obs[batch].copy()
    
    # set 'cell_type' and 'batch' for scGen
    adata.obs['cell_type'] = adata.obs[cell_type].copy()
    adata.obs['batch'] = adata.obs[batch].copy()
    
    # feature selection
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    adata = adata[:, adata.var['highly_variable'] == True].copy()
    
    network = scgen.VAEArith(x_dimension= adata.shape[1], model_path=model_path)
    network.train(train_data=adata, n_epochs=epochs)
    corrected_adata = scgen.batch_removal(network, adata)
    network.sess.close()
    
    # reset fields (just in case they were overwritten)
    adata.obs[cell_type] = cell_types
    adata.obs[batch] = batches
    
    return corrected_adata

def runSeurat(adata, batch, hvg=None):
    checkSanity(adata, batch, hvg)
    ro.r('library(Seurat)')
    ro.r('library(scater)')
    anndata2ri.activate()
    
    tmp = anndata.AnnData(X=adata.X.sorted_indices(), obs=adata.obs)
    ro.globalenv['adata'] = tmp
    ro.r('sobj = as.Seurat(adata, counts=NULL, data = "X")')

    ro.r(f'batch_list = SplitObject(sobj, split.by = "{batch}")')
    #ro.r('to_integrate <- Reduce(intersect, lapply(batch_list, rownames))')
    ro.r('anchors = FindIntegrationAnchors('+
        'object.list = batch_list, '+
        'anchor.features = 2000,'+
        'scale = T,'+
        'l2.norm = T,'+
        'dims = 1:30,'+
        'k.anchor = 5,'+
        'k.filter = 200,'+
        'k.score = 30,'+
        'max.features = 200,'+
        'eps = 0)'
    )
    ro.r('integrated = IntegrateData('+
        'anchorset = anchors,'+
        'new.assay.name = "integrated",'+
        'features = NULL,'+
        'features.to.integrate = NULL,'+
        'dims = 1:30,'+
        'k.weight = 100,'+
        'weight.reduction = NULL,'+
        'sd.weight = 1,'+
        'sample.tree = NULL,'+
        'preserve.order = F,'+
        'do.cpp = T,'+
        'eps = 0,'+
        'verbose = T)'
    )
    integrated = ro.r('as.SingleCellExperiment(integrated)')
    anndata2ri.deactivate()
    return integrated

def runHarmony(adata, batch, hvg = None):
    checkSanity(adata, batch, hvg)
    #import_rpy2()
    ro.pandas2ri.activate()
    ro.r('library(harmony)')

    pca = sc.pp.pca(adata, svd_solver='arpack', copy=True).obsm['X_pca']
    method = adata.obs[batch]

    ro.globalenv['pca'] = pca
    ro.globalenv['method'] = method

    ro.r(f'harmonyEmb <- HarmonyMatrix(pca, method, "{batch}", do_pca= F)')
    emb = ro.r('harmonyEmb')
    ro.pandas2ri.deactivate()
    out = adata.copy()
    out.obsm['X_pca']= emb
    out.uns['emb']=True

    return out

def runMNN(adata, batch, hvg = None):
    import mnnpy
    checkSanity(adata, batch, hvg)
    split = splitBatches(adata, batch)

    corrected = mnnpy.mnn_correct(*split, var_subset=hvg)

    return corrected[0]

def runBBKNN(adata, batch, hvg=None):
    import bbknn
    checkSanity(adata, batch, hvg)
    sc.pp.pca(adata, svd_solver='arpack')
    corrected = bbknn.bbknn(adata, batch_key=batch, copy=True)
    return corrected

def runConos(adata, batch, hvg=None):
    checkSanity(adata, batch, hvg)
    anndata2ri.activate()
    ro.r('library(Seurat)')
    ro.r('library(scater)')
    ro.r('library(conos)')

    ro.globalenv['adata_c'] = adata
    ro.r('sobj = as.Seurat(adata_c, counts = "counts", data = "X")')
    ro.r(f'batch_list = SplitObject(sobj, split.by = "{batch}")')

    ro.r('con <- Conos$new(batch_list)')
    ro.r('con$buildGraph(k=15, k.self=5, space="PCA", ncomps=30)')
    os.mkdir('conos_tmp')
    ro.r('saveConosForScanPy(con, output.path="conos_tmp/", verbose=T))')

    DATA_PATH = os.path.expanduser('conos_tmp/')
    pca_df = pd.read_csv(DATA_PATH+'pca.csv')
    
    graph_conn_mtx = sp.io.mmread(DATA_PATH + "graph_connectivities.mtx")
    graph_dist_mtx = sp.io.mmread(DATA_PATH + "graph_distances.mtx")
    
    out = adata.copy()
    
    out.X_pca = pca_df.values
    
    out.uns['neighbors'] = dict(connectivities=graph_conn_mtx.tocsr(), distances=graph_dist_mtx.tocsr())
    
    anndata2ri.deactivate()
    return out

    


if __name__=="__main__":
    adata = sc.read('testing.h5ad')
    #emb, corrected = runScanorama(adata, 'method', False)
    #print(emb)
    #print(corrected)


        

