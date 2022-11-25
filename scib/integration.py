import os
import tempfile

import numpy as np
import scanpy as sc
import scipy as sp
from scipy.sparse import issparse

from . import utils
from .exceptions import OptionalDependencyNotInstalled


def harmony(adata, batch, hvg=None, **kwargs):
    """Harmony wrapper function

    Based on `harmony-pytorch <https://github.com/lilab-bcb/harmony-pytorch>`_ version 0.1.7

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix as well as an embedding representation of the
        corrected data
    """
    try:
        from harmony import harmonize
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

    utils.check_sanity(adata, batch, hvg)
    sc.tl.pca(adata)
    adata.obsm["X_emb"] = harmonize(adata.obsm["X_pca"], adata.obs, batch_key=batch)

    return adata


def scanorama(adata, batch, hvg=None, **kwargs):
    """Scanorama wrapper function

    Based on `scanorama <https://github.com/brianhie/scanorama>`_ version 1.7.0

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix as well as an embedding representation of the
        corrected data
    """
    try:
        import scanorama
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

    utils.check_sanity(adata, batch, hvg)
    split, categories = utils.split_batches(adata.copy(), batch, return_categories=True)
    corrected = scanorama.correct_scanpy(split, return_dimred=True, **kwargs)
    corrected = utils.merge_adata(
        *corrected, batch_key=batch, batch_categories=categories, index_unique=None
    )
    corrected.obsm["X_emb"] = corrected.obsm["X_scanorama"]
    # corrected.uns['emb']=True

    return corrected


def trvae(adata, batch, hvg=None):
    """trVAE wrapper function

    Based on `trVAE <https://github.com/theislab/trVAE>`_ version 1.1.2

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix as well as an embedding representation of the
        corrected data
    """
    try:
        import trvae
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

    utils.check_sanity(adata, batch, hvg)
    n_batches = len(adata.obs[batch].cat.categories)

    train_adata, valid_adata = trvae.utils.train_test_split(adata, train_frac=0.80)

    condition_encoder = trvae.utils.create_dictionary(
        adata.obs[batch].cat.categories.tolist(), []
    )

    network = trvae.archs.trVAEMulti(
        x_dimension=train_adata.shape[1],
        n_conditions=n_batches,
        output_activation="relu",
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

    adata.obsm["X_emb"] = adata.obsm["mmd_latent"]
    del adata.obsm["mmd_latent"]
    adata.X = adata.obsm["reconstructed"]

    return adata


def trvaep(adata, batch, hvg=None):
    """trVAE wrapper function (``pytorch`` implementatioon)

    Based on `trvaep`_ version 0.1.0

    .. _trvaep: https://github.com/theislab/trvaep

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix as well as an embedding representation of the
        corrected data
    """
    try:
        import trvaep
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

    utils.check_sanity(adata, batch, hvg)
    n_batches = adata.obs[batch].nunique()

    # Densify the data matrix
    if issparse(adata.X):
        adata.X = adata.X.A

    model = trvaep.CVAE(
        adata.n_vars,
        num_classes=n_batches,
        encoder_layer_sizes=[64, 32],
        decoder_layer_sizes=[32, 64],
        latent_dim=10,
        alpha=0.0001,
        use_mmd=True,
        beta=1,
        output_activation="ReLU",
    )

    # Note: set seed for reproducibility of results
    trainer = trvaep.Trainer(model, adata, condition_key=batch, seed=42)

    trainer.train_trvae(300, 1024, early_patience=50)

    # Get the dominant batch covariate
    main_batch = adata.obs[batch].value_counts().idxmax()

    # Get latent representation
    latent_y = model.get_y(
        adata.X,
        c=model.label_encoder.transform(np.tile(np.array([main_batch]), len(adata))),
    )
    adata.obsm["X_emb"] = latent_y

    # Get reconstructed feature space:
    data = model.predict(x=adata.X, y=adata.obs[batch].tolist(), target=main_batch)
    adata.X = data

    return adata


def scgen(adata, batch, cell_type, epochs=100, hvg=None, **kwargs):
    """scGen wrapper function

    Based on `scgen`_ version 2.1.0 with parametrization taken from the tutorial `notebook`_.

    .. _scgen: https://github.com/theislab/scgen
    .. _notebook: https://scgen.readthedocs.io/en/stable/tutorials/scgen_batch_removal.html

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix
    """
    try:
        from scgen import SCGEN
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

    utils.check_sanity(adata, batch, hvg)

    net_adata = adata.copy()
    if hvg is not None:
        net_adata = net_adata[:, hvg].copy()

    SCGEN.setup_anndata(net_adata, batch_key=batch, labels_key=cell_type)
    model = SCGEN(net_adata)
    model.train(
        max_epochs=epochs,
        batch_size=32,
        early_stopping=True,
        early_stopping_patience=25,
    )
    corrected_adata = model.batch_removal(**kwargs)
    return corrected_adata


def scvi(adata, batch, hvg=None, return_model=False, max_epochs=None):
    """scVI wrapper function

    Based on scvi-tools version >=0.16.0 (available through `conda <https://docs.scvi-tools.org/en/stable/installation.html>`_)

    .. note::
        scVI expects only non-normalized (count) data on highly variable genes!

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix as well as an embedding representation of the
        corrected data
    """
    try:
        from scvi.model import SCVI
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

    utils.check_sanity(adata, batch, hvg)

    # Check for counts data layer
    if "counts" not in adata.layers:
        raise TypeError(
            "Adata does not contain a `counts` layer in `adata.layers[`counts`]`"
        )

    # Defaults from SCVI github tutorials scanpy_pbmc3k and harmonization
    n_latent = 30
    n_hidden = 128
    n_layers = 2

    # copying to not return values added to adata during setup_anndata
    net_adata = adata.copy()
    if hvg is not None:
        net_adata = adata[:, hvg].copy()
    SCVI.setup_anndata(net_adata, layer="counts", batch_key=batch)

    vae = SCVI(
        net_adata,
        gene_likelihood="nb",
        n_layers=n_layers,
        n_latent=n_latent,
        n_hidden=n_hidden,
    )
    train_kwargs = {"train_size": 1.0}
    if max_epochs is not None:
        train_kwargs["max_epochs"] = max_epochs
    vae.train(**train_kwargs)
    adata.obsm["X_emb"] = vae.get_latent_representation()

    if not return_model:
        return adata
    else:
        return vae


def scanvi(adata, batch, labels, hvg=None, max_epochs=None):
    """scANVI wrapper function

    Based on scvi-tools version >=0.16.0 (available through `conda <https://docs.scvi-tools.org/en/stable/installation.html>`_)

    .. note::
        Use non-normalized (count) data for scANVI!

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param labels: label key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix as well as an embedding representation of the
        corrected data
    """
    try:
        from scvi.model import SCANVI
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

    # # Defaults from SCVI github tutorials scanpy_pbmc3k and harmonization
    # this n_epochs_scVI is now default in scvi-tools
    if max_epochs is None:
        n_epochs_scVI = int(np.min([round((20000 / adata.n_obs) * 400), 400]))  # 400
        n_epochs_scANVI = int(np.min([10, np.max([2, round(n_epochs_scVI / 3.0)])]))
    else:
        n_epochs_scVI = max_epochs
        n_epochs_scANVI = max_epochs

    vae = scvi(adata, batch, hvg, return_model=True, max_epochs=n_epochs_scVI)

    # STEP 2: RUN scVI to initialize scANVI
    scanvae = SCANVI.from_scvi_model(
        scvi_model=vae,
        labels_key=labels,
        unlabeled_category="UnknownUnknown",  # pick anything definitely not in a dataset
    )
    scanvae.train(max_epochs=n_epochs_scANVI, train_size=1.0)
    adata.obsm["X_emb"] = scanvae.get_latent_representation()

    return adata


def mnn(adata, batch, hvg=None, **kwargs):
    """MNN wrapper function (``mnnpy`` implementation)

    Based on `mnnpy package <https://github.com/chriscainx/mnnpy>`_ version 0.1.9.5

    .. note:

        ``mnnpy`` might break with newer versions of ``numpy`` and ``pandas``

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix
    """
    try:
        import mnnpy
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

    utils.check_sanity(adata, batch, hvg)
    split, categories = utils.split_batches(adata, batch, return_categories=True)

    corrected, _, _ = mnnpy.mnn_correct(
        *split,
        var_subset=hvg,
        batch_key=batch,
        batch_categories=categories,
        index_unique=None,
        **kwargs,
    )

    return corrected


def bbknn(adata, batch, hvg=None, **kwargs):
    """BBKNN wrapper function

    Based on `bbknn package <https://github.com/Teichlab/bbknn>`_ version 1.3.9

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :params \\**kwargs: additional parameters for BBKNN
    :return: ``anndata`` object containing the corrected graph
    """
    try:
        import bbknn
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

    utils.check_sanity(adata, batch, hvg)
    sc.pp.pca(adata, svd_solver="arpack")
    if adata.n_obs < 1e5:
        return bbknn.bbknn(adata, batch_key=batch, copy=True, **kwargs)
    if adata.n_obs >= 1e5:
        return bbknn.bbknn(
            adata, batch_key=batch, neighbors_within_batch=25, copy=True, **kwargs
        )


def saucie(adata, batch):
    """SAUCIE wrapper function

    Using SAUCIE `source code <https://github.com/KrishnaswamyLab/SAUCIE>`_.
    Parametrisation from https://github.com/KrishnaswamyLab/SAUCIE/blob/master/scripts/SAUCIE.py

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :return: ``anndata`` object containing the corrected embedding
    """
    try:
        import SAUCIE
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

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
    ret.obsm["X_emb"] = saucie.get_reconstruction(loader_eval)[0]
    ret.X = pca_op.inverse_transform(ret.obsm["X_emb"])

    return ret


def combat(adata, batch):
    """ComBat wrapper function (``scanpy`` implementation)

    Using scanpy implementation of `Combat <https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.combat.html>`_

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :return: ``anndata`` object containing the corrected feature matrix
    """
    adata_int = adata.copy()
    sc.pp.combat(adata_int, key=batch)
    return adata_int


def desc(adata, batch, res=0.8, ncores=None, tmp_dir=None, use_gpu=False):
    """DESC wrapper function

    Based on `desc package <https://github.com/eleozzr/desc>`_ version 2.0.3.
    Parametrization was taken from: https://github.com/eleozzr/desc/issues/28 as suggested by the developer (rather
    than from the tutorial notebook).

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :return: ``anndata`` object containing the corrected embedding
    """
    try:
        import desc
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

    if tmp_dir is None:
        temp_dir = tempfile.TemporaryDirectory()
        tmp_dir = temp_dir.name

    # Set number of CPUs to all available
    if ncores is None:
        ncores = os.cpu_count()

    adata_out = adata.copy()

    adata_out = desc.scale_bygroup(adata_out, groupby=batch, max_value=6)

    adata_out = desc.train(
        adata_out,
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
        do_umap=False,
    )

    adata_out.obsm["X_emb"] = adata_out.obsm["X_Embeded_z" + str(res)]

    return adata_out
