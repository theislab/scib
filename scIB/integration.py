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
from scipy.sparse import issparse

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


def runTrVaep(adata, batch, hvg=None):
    checkSanity(adata, batch, hvg)
    import trvaep

    n_batches = adata.obs[batch].nunique()
    
    # Densify the data matrix
    if issparse(adata.X):
        adata.X = adata.X.A

    model = trvaep.CVAE(adata.n_vars, num_classes=n_batches,
                        encoder_layer_sizes=[64, 32],
                        decoder_layer_sizes=[32, 64], latent_dim=10,
                        alpha=0.0001, use_mmd=True, beta=1,
                        output_activation="ReLU")
    
    # Note: set seed for reproducibility of results
    trainer = trvaep.Trainer(model, adata, condition_key=batch, seed=42)
    
    trainer.train_trvae(300, 1024, early_patience=50)

    # Get the dominant batch covariate
    main_batch = adata.obs[batch].value_counts().idxmax()

    # Get latent representation
    latent_y = model.get_y(
        adata.X, c=model.label_encoder.transform(
            np.tile(np.array([main_batch]), len(adata))))
    adata.obsm['X_emb'] = latent_y

    # Get reconstructed feature space:
    data = model.predict(x=adata.X, y=adata.obs[batch].tolist(),
                         target=main_batch)
    adata.X = data

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
    n_epochs=np.min([round((20000/adata.n_obs)*400), 400])
    n_latent=30
    n_hidden=128
    n_layers=2
    
    net_adata = adata.copy()
    net_adata.X = adata.layers['counts']
    del net_adata.layers['counts']
    # Ensure that the raw counts are not accidentally used
    del net_adata.raw # Note that this only works from anndata 0.7

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
        train_size=1.0,
        use_cuda=False,
    )

    trainer.train(n_epochs=n_epochs, lr=1e-3)

    full = trainer.create_posterior(trainer.model, net_adata, indices=np.arange(len(net_adata)))
    latent, _, _ = full.sequential().get_latent()

    adata.obsm['X_emb'] = latent

    return adata

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
    if adata.n_obs <1e5:
        corrected = bbknn.bbknn(adata, batch_key=batch, copy=True)
    if adata.n_obs >=250000:
        corrected = bbknn.bbknn(adata, batch_key=batch, neighbors_within_batch=25, copy=True)
    return corrected


def runSaucie(adata, batch):
    """
    parametrisation from https://github.com/KrishnaswamyLab/SAUCIE/blob/master/scripts/SAUCIE.py
    """
    import SAUCIE
    import sklearn.decomposition
    pca_op = sklearn.decomposition.PCA(100)
    expr = adata.X.todense()
    data = pca_op.fit_transform(expr)
    saucie = SAUCIE.SAUCIE(100, lambda_b=0.1)
    loader_train = SAUCIE.Loader(data, labels=adata.obs[batch].cat.codes, shuffle=True)
    loader_eval = SAUCIE.Loader(data, labels=adata.obs[batch].cat.codes, shuffle=False)
    saucie.train(loader_train, steps=5000)
    ret = adata.copy()
    ret.obsm['X_emb'] = saucie.get_reconstruction(loader_eval)[0]
    ret.X = pca_op.inverse_transform(ret.obsm['X_emb'])
                                                    
    return ret
                                                             
    

def runCombat(adata, batch):
    adata_int = adata.copy()
    sc.pp.combat(adata_int, key=batch)
    return adata_int


def runDESC(adata, batch, res=0.8, ncores=None, tmp_dir='/localscratch/'):
    """
    Convenience function to run DESC. Parametrization was taken from:
    https://github.com/eleozzr/desc/issues/28
    as suggested by the developer (rather than from the tutorial notebook).
    """
    import desc

    # Set number of CPUs to all available
    if ncores is None:
        ncores = os.cpu_count()

    adata_out = adata.copy()

    adata_out = desc.scale_bygroup(adata_out, groupby=batch, max_value=6)
    
    adata_out = desc.train(adata_out,
                     dims=[adata.shape[1],128,32],
                     tol=0.001,
                     n_neighbors=10,
                     batch_size=256,
                     louvain_resolution=res,
                     save_dir=tmp_dir,
                     do_tsne=False,
                     use_GPU=False,
                     num_Cores=ncores,
                     save_encoder_weights=False,
                     use_ae_weights=False,
                     do_umap=False)
    
    adata.obsm['X_emb'] = adata_out.obsm['X_Embeded_z'+str(res)]

    return adata


if __name__=="__main__":
    adata = sc.read('testing.h5ad')
    #emb, corrected = runScanorama(adata, 'method', False)
    #print(emb)
    #print(corrected)


        

