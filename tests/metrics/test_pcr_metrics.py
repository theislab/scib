from tests.common import *
import pandas as pd


def test_pc_regression(adata_pbmc):
    scIB.me.pc_regression(adata_pbmc.X, adata_pbmc.obs["batch"])


def test_pcr_batch(adata_pbmc, verbose=True):
    # no PCA precomputed
    adata_int = adata_pbmc.copy()

    score = scIB.me.pcr_comparison(
        adata_pbmc,
        adata_int,
        variable='batch',
        n_comps=50,
        scale=True,
        verbose=verbose
    )
    LOGGER.info(f"no PCA precomputed: {score}")
    assert score == 0


def test_pcr_batch_precomputed(adata, adata_pca):
    adata = adata_pca(adata)
    adata_int = adata.copy()
    score = scIB.me.pcr_comparison(adata, adata_int, variable='batch', scale=True)
    LOGGER.info(f"precomputed PCA: {score}")
    assert score == 0  # same PCA values -> difference should be 0


def test_pcr_batch_embedding(adata_pbmc, embed_factory):
    # use different embedding
    adata_int = embed_factory(adata_pbmc, type_='full')
    score = scIB.me.pcr_comparison(
        adata_pbmc, adata_int,
        variable='batch',
        embed='X_emb',
        n_comps=50,
        scale=True
    )
    LOGGER.info(f"using embedding: {score}")
    assert score >= 0
    assert score <= 1
    assert score < 1e-6
