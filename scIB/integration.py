#!/bin/env python

### D. C. Strobl, M. Müller; 2019-07-23

""" This module provides a toolkit for running a large range of single cell data integration methods
    as well as tools and metrics to benchmark them.
"""

import scipy as sp
from scIB.utils import *
import os
import anndata
from scIB.exceptions import IntegrationMethodNotFound

import rpy2.rinterface_lib.callbacks
import logging

rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)  # Ignore R warning messages
from scipy.sparse import issparse


# functions for running the methods

def runScanorama(adata, batch, hvg=None):
    try:
        import scanorama
    except ModuleNotFoundError as e:
        raise IntegrationMethodNotFound(e)

    checkSanity(adata, batch, hvg)
    split, categories = splitBatches(adata.copy(), batch, return_categories=True)
    corrected = scanorama.correct_scanpy(split, return_dimred=True)
    corrected = anndata.AnnData.concatenate(
        *corrected, batch_key=batch, batch_categories=categories, index_unique=None
    )
    corrected.obsm['X_emb'] = corrected.obsm['X_scanorama']
    # corrected.uns['emb']=True

    return corrected


def runTrVae(adata, batch, hvg=None):
    try:
        import trvae
    except ModuleNotFoundError as e:
        raise IntegrationMethodNotFound(e)

    checkSanity(adata, batch, hvg)
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
    try:
        import trvaep
    except ModuleNotFoundError as e:
        raise IntegrationMethodNotFound(e)

    checkSanity(adata, batch, hvg)
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


def runScGen(adata, batch, cell_type, epochs=100, hvg=None, model_path='/localscratch'):
    """
    Parametrization taken from the tutorial notebook at:
    https://nbviewer.jupyter.org/github/M0hammadL/scGen_notebooks/blob/master/notebooks/scgen_batch_removal.ipynb
    """
    try:
        import scgen
    except ModuleNotFoundError as e:
        raise IntegrationMethodNotFound(e)

    checkSanity(adata, batch, hvg)

    # Fit the model
    network = scgen.VAEArith(x_dimension=adata.shape[1], model_path=model_path)
    network.train(train_data=adata, n_epochs=epochs, save=False)
    corrected_adata = scgen.batch_removal(network, adata, batch_key=batch, cell_label_key=cell_type)

    network.sess.close()

    return corrected_adata


def runScvi(adata, batch, hvg=None):
    # Use non-normalized (count) data for scvi!
    # Expects data only on HVGs
    try:
        from scvi.models import VAE
        from scvi.inference import UnsupervisedTrainer
        from sklearn.preprocessing import LabelEncoder
        from scvi.dataset import AnnDatasetFromAnnData
    except ModuleNotFoundError as e:
        raise IntegrationMethodNotFound(e)

    checkSanity(adata, batch, hvg)

    # Check for counts data layer
    if 'counts' not in adata.layers:
        raise TypeError('Adata does not contain a `counts` layer in `adata.layers[`counts`]`')

    # Defaults from SCVI github tutorials scanpy_pbmc3k and harmonization
    n_epochs = np.min([round((20000 / adata.n_obs) * 400), 400])
    n_latent = 30
    n_hidden = 128
    n_layers = 2

    net_adata = adata.copy()
    net_adata.X = adata.layers['counts']
    del net_adata.layers['counts']
    # Ensure that the raw counts are not accidentally used
    del net_adata.raw  # Note that this only works from anndata 0.7

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


def runScanvi(adata, batch, labels):
    # Use non-normalized (count) data for scanvi!
    try:
        from scvi.models import VAE, SCANVI
        from scvi.inference import UnsupervisedTrainer, SemiSupervisedTrainer
        from sklearn.preprocessing import LabelEncoder
        from scvi.dataset import AnnDatasetFromAnnData
    except ModuleNotFoundError as e:
        raise IntegrationMethodNotFound(e)

    import numpy as np

    if 'counts' not in adata.layers:
        raise TypeError('Adata does not contain a `counts` layer in `adata.layers[`counts`]`')
    # Check for counts data layer

    # STEP 1: prepare the data
    net_adata = adata.copy()
    net_adata.X = adata.layers['counts']
    del net_adata.layers['counts']
    # Ensure that the raw counts are not accidentally used
    del net_adata.raw  # Note that this only works from anndata 0.7

    # Define batch indices
    le = LabelEncoder()
    net_adata.obs['batch_indices'] = le.fit_transform(net_adata.obs[batch].values)
    net_adata.obs['labels'] = le.fit_transform(net_adata.obs[labels].values)

    net_adata = AnnDatasetFromAnnData(net_adata)

    print("scANVI dataset object with {} batches and {} cell types".format(net_adata.n_batches, net_adata.n_labels))

    # if hvg is True:
    #    # this also corrects for different batches by default
    #    net_adata.subsample_genes(2000, mode="seurat_v3")

    # # Defaults from SCVI github tutorials scanpy_pbmc3k and harmonization
    n_epochs_scVI = np.min([round((20000 / adata.n_obs) * 400), 400])  # 400
    n_epochs_scANVI = int(np.min([10, np.max([2, round(n_epochs_scVI / 3.)])]))
    n_latent = 30
    n_hidden = 128
    n_layers = 2

    # STEP 2: RUN scVI to initialize scANVI

    vae = VAE(
        net_adata.nb_genes,
        reconstruction_loss='nb',
        n_batch=net_adata.n_batches,
        n_latent=n_latent,
        n_hidden=n_hidden,
        n_layers=n_layers,
    )

    trainer = UnsupervisedTrainer(
        vae,
        net_adata,
        train_size=1.0,
        use_cuda=False,
    )

    trainer.train(n_epochs=n_epochs_scVI, lr=1e-3)

    # STEP 3: RUN scANVI

    scanvi = SCANVI(net_adata.nb_genes, net_adata.n_batches, net_adata.n_labels,
                    n_hidden=n_hidden, n_latent=n_latent, n_layers=n_layers, dispersion='gene',
                    reconstruction_loss='nb')
    scanvi.load_state_dict(trainer.model.state_dict(), strict=False)

    # use default parameter from semi-supervised trainer class
    trainer_scanvi = SemiSupervisedTrainer(scanvi, net_adata)
    # use all cells as labelled set
    trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(trainer_scanvi.model, net_adata,
                                                                  indices=np.arange(len(net_adata)))
    # put one cell in the unlabelled set
    trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=[0])
    trainer_scanvi.train(n_epochs=n_epochs_scANVI)

    # extract info from posterior
    scanvi_full = trainer_scanvi.create_posterior(trainer_scanvi.model, net_adata, indices=np.arange(len(net_adata)))
    latent, _, _ = scanvi_full.sequential().get_latent()

    adata.obsm['X_emb'] = latent

    return adata


def runMNN(adata, batch, hvg=None):
    try:
        import mnnpy
    except ModuleNotFoundError as e:
        raise IntegrationMethodNotFound(e)

    checkSanity(adata, batch, hvg)
    split, categories = splitBatches(adata, batch, return_categories=True)

    corrected, _, _ = mnnpy.mnn_correct(
        *split, var_subset=hvg, batch_key=batch, batch_categories=categories, index_unique=None
    )

    return corrected


def runBBKNN(adata, batch, hvg=None):
    try:
        import bbknn
    except ModuleNotFoundError as e:
        raise IntegrationMethodNotFound(e)

    checkSanity(adata, batch, hvg)
    sc.pp.pca(adata, svd_solver='arpack')
    if adata.n_obs < 1e5:
        return bbknn.bbknn(adata, batch_key=batch, copy=True)
    if adata.n_obs >= 1e5:
        return bbknn.bbknn(adata, batch_key=batch, neighbors_within_batch=25, copy=True)


def runSaucie(adata, batch):
    """
    parametrisation from https://github.com/KrishnaswamyLab/SAUCIE/blob/master/scripts/SAUCIE.py
    """
    try:
        import SAUCIE
    except ModuleNotFoundError as e:
        raise IntegrationMethodNotFound(e)

    import sklearn.decomposition
    pca_op = sklearn.decomposition.PCA(100)
    if isinstance(adata.X, sp.sparse.csr_matrix):
        expr = adata.X.A
    else:
        expr = adata.X
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


def runDESC(adata, batch, res=0.8, ncores=None, tmp_dir='/localscratch/tmp_desc/', use_gpu=False):
    """
    Convenience function to run DESC. Parametrization was taken from:
    https://github.com/eleozzr/desc/issues/28
    as suggested by the developer (rather than from the tutorial notebook).
    """
    try:
        import desc
    except ModuleNotFoundError as e:
        raise IntegrationMethodNotFound(e)

    # Set number of CPUs to all available
    if ncores is None:
        ncores = os.cpu_count()

    adata_out = adata.copy()

    adata_out = desc.scale_bygroup(adata_out, groupby=batch, max_value=6)

    adata_out = desc.train(adata_out,
                           dims=[adata.shape[1], 128, 32],
                           tol=0.001,
                           n_neighbors=10,
                           batch_size=256,
                           louvain_resolution=res,
                           save_encoder_weights=False,
                           save_dir=tmp_dir,
                           do_tsne=False,
                           use_GPU=use_gpu,
                           num_Cores=ncores,
                           use_ae_weights=False,
                           do_umap=False)

    adata_out.obsm['X_emb'] = adata_out.obsm['X_Embeded_z' + str(res)]

    return adata_out
