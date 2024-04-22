from scipy.sparse import csr_matrix

import scib
from tests.common import LOGGER, add_embed, assert_near_exact


def test_pc_regression(adata):
    score = scib.me.pc_regression(adata.X, adata.obs["batch"])
    LOGGER.info(score)
    assert_near_exact(score, 0, diff=1e-4)


def test_pc_regression_sparse(adata):
    score = scib.me.pc_regression(
        csr_matrix(adata.X),
        covariate=adata.obs["batch"],
        n_comps=adata.n_vars,
    )
    LOGGER.info(score)
    assert_near_exact(score, 0, diff=1e-4)


def test_pcr_sklearn(adata_pca):
    score = scib.me.pcr(
        adata_pca, covariate="celltype", linreg_method="sklearn", verbose=True
    )
    LOGGER.info(score)
    assert_near_exact(score, 0.3371261556141021, diff=1e-3)


def test_pcr_numpy(adata_pca):
    score = scib.me.pcr(
        adata_pca, covariate="celltype", linreg_method="numpy", verbose=True
    )
    LOGGER.info(score)
    assert_near_exact(score, 0.3371261556141021, diff=1e-3)


def test_pcr_implementations(adata_pca):
    score_sklearn = scib.me.pcr(
        adata_pca,
        covariate="celltype",
        linreg_method="sklearn",
    )
    LOGGER.info(f"sklearn score: {score_sklearn}")

    score_numpy = scib.me.pcr(
        adata_pca,
        covariate="celltype",
        linreg_method="numpy",
    )
    LOGGER.info(f"numpy score: {score_numpy}")

    assert_near_exact(score_sklearn, score_numpy, diff=1e-3)


def test_pcr_comparison_batch(adata):
    score = scib.me.pcr_comparison(
        adata, adata, covariate="batch", n_comps=50, scale=True
    )
    LOGGER.info(f"no PCA precomputed: {score}")
    assert_near_exact(score, 0, diff=1e-6)


def test_pcr_comparison_batch_precomputed(adata_pca):
    score = scib.me.pcr_comparison(adata_pca, adata_pca, covariate="batch", scale=True)
    LOGGER.info(f"precomputed PCA: {score}")
    assert_near_exact(score, 0, diff=1e-6)


def test_pcr_comparison_batch_embedding(adata):
    # use different embedding
    score = scib.me.pcr_comparison(
        adata_pre=adata,
        adata_post=add_embed(adata, type_="full"),
        covariate="batch",
        embed="X_emb",
        n_comps=50,
        scale=True,
    )
    LOGGER.info(f"using X_emb: {score}")
    assert_near_exact(score, 0, diff=1e-6)
