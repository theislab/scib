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
    corrected.obsm['X_emb']= emb
    #corrected.uns['emb']=True

    return corrected

def runTrVae(adata, batch, hvg=None):
    checkSanity(adata, batch, hvg)
    import trvae

    n_batches = len(adata.obs[batch].cat.categories)

    train_adata, valid_adata = trvae.utils.train_test_split(
        adata,
        train_frac=0.80
    )

    condition_encoder = trvae.utils.create_dictionary(
        adata.obs[batch].cat.categories.tolist(), [])

    network = trvae.archs.trVAEMulti(
        x_dimension=train_adata.shape[1],
        n_conditions=n_batches,
        output_activation='relu'
    )

    network.train(
        train_adata,
        valid_adata,
        condition_key=batch,
        condition_encoder=condition_encoder,
        verbose=0,
    )

    labels, _ = trvae.tl.label_encoder(
        adata,
        condition_key=batch,
        label_encoder=condition_encoder,
    )

    network.get_corrected(adata, labels, return_z=False)
    
    adata.obsm['X_emb'] = adata.obsm['mmd_latent']
    del adata.obsm['mmd_latent']
    adata.X = adata.obsm['reconstructed']
    
    return adata


def runScvi(adata, batch, hvg=None):
    # Use non-normalized (count) data for scvi!
    # Expects data only on HVGs
    
    checkSanity(adata, batch, hvg)

    # Check for counts data layer
    if 'counts' not in adata.layers:
        raise TypeError('Adata does not contain a `counts` layer in `adata.layers[`counts`]`')

    from scvi.models import VAE
    from scvi.inference import UnsupervisedTrainer
    from sklearn.preprocessing import LabelEncoder
    from scvi.dataset import AnnDatasetFromAnnData

    # Defaults from SCVI github tutorials scanpy_pbmc3k and harmonization
    n_epochs=100
    n_latent=30
    n_hidden=128
    n_layers=2
    
    net_adata = adata.copy()
    net_adata.X = adata.layers['counts']
    del net_adata.layers['counts']
    net_adata.raw = None # Ensure that the raw counts are not accidentally used

    # Define batch indices
    le = LabelEncoder()
    net_adata.obs['batch_indices'] = le.fit_transform(net_adata.obs[batch].values)

    net_adata = AnnDatasetFromAnnData(net_adata)

    vae = VAE(
        net_adata.nb_genes,
        reconstruction_loss='nb',
        n_batch=net_adata.n_batches,
        n_layers=n_layers,
        n_latent=n_latent,
        n_hidden=n_hidden,
    )

    trainer = UnsupervisedTrainer(
        vae,
        net_adata,
        train_size=1,
        use_cuda=False,
    )

    trainer.train(n_epochs=n_epochs, lr=1e-3)

    full = trainer.create_posterior(trainer.model, net_adata, indices=np.arange(len(net_adata)))
    latent, _, _ = full.sequential().get_latent()

    adata.obsm['X_emb'] = latent

    return adata


def runSeurat(adata, batch, hvg=None):
    checkSanity(adata, batch, hvg)
    ro.r('library(Seurat)')
    ro.r('library(scater)')
    anndata2ri.activate()
    
    #tmp = anndata.AnnData(X=adata.X.sorted_indices(), obs=adata.obs)
    #ro.globalenv['adata'] = tmp

    ro.globalenv['adata'] = adata
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
    out.obsm['X_emb']= emb
    #out.uns['emb']=True

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


        

