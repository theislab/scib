from scipy.sparse import csr_matrix

import scib
from tests.common import LOGGER, add_embed, assert_near_exact


def test_pc_regression(adata):
    score = scib.me.pc_regression(adata.X, adata.obs["batch"])
    LOGGER.info(f"using embedding: {score}")
    assert_near_exact(score, 0, diff=1e-4)


def test_pc_regression_sparse(adata):
    x = csr_matrix(adata.X)
    score = scib.me.pc_regression(x, adata.obs["batch"], n_comps=x.shape[1])
    LOGGER.info(f"using embedding: {score}")
    assert_near_exact(score, 0, diff=1e-4)


def test_pcr_batch(adata):
    # no PCA precomputed
    score = scib.me.pcr_comparison(
        adata, adata, covariate="batch", n_comps=50, scale=True
    )
    LOGGER.info(f"no PCA precomputed: {score}")
    assert_near_exact(score, 0, diff=1e-6)


def test_pcr_batch_precomputed(adata_pca):
    score = scib.me.pcr_comparison(adata_pca, adata_pca, covariate="batch", scale=True)
    LOGGER.info(f"precomputed PCA: {score}")
    assert_near_exact(score, 0, diff=1e-6)


def test_pcr_batch_embedding(adata):
    # use different embedding
    score = scib.me.pcr_comparison(
        adata_pre=adata,
        adata_post=add_embed(adata, type_="full"),
        covariate="batch",
        embed="X_emb",
        n_comps=50,
        scale=True,
    )
    LOGGER.info(f"using embedding: {score}")
    assert_near_exact(score, 0, diff=1e-6)
