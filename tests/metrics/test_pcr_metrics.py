import pytest
from scipy.sparse import csr_matrix

import scib
from tests.common import LOGGER, add_embed, assert_near_exact


@pytest.mark.parametrize("sparse", [False, True])
def test_pc_regression(adata, sparse):
    if sparse:
        adata.X = csr_matrix(adata.X)
    score = scib.me.pc_regression(
        adata.X,
        covariate=adata.obs["batch"],
        n_comps=adata.n_vars,
    )
    LOGGER.info(score)
    assert_near_exact(score, 0, diff=1e-3)


@pytest.mark.parametrize("linreg_method", ["numpy", "sklearn"])
@pytest.mark.parametrize("n_threads", [1, 2])
def test_pcr_timing(adata_pca, n_threads, linreg_method):
    import timeit

    import anndata as ad
    import scanpy as sc

    # scale up anndata
    adata = ad.concat([adata_pca] * 100)
    print(f"compute PCA on {adata.n_obs} cells...")
    sc.pp.pca(adata)

    runs = 10
    timing = timeit.timeit(
        lambda: scib.me.pcr(
            adata,
            covariate="celltype",
            linreg_method=linreg_method,
            verbose=False,
            n_threads=n_threads,
        ),
        number=runs,
    )
    LOGGER.info(f"timeit: {timing:.2f}s for {runs} runs")

    # test pcr value
    score = scib.me.pcr(
        adata_pca,
        covariate="celltype",
        linreg_method=linreg_method,
        verbose=True,
        n_threads=1,
    )
    LOGGER.info(score)
    assert_near_exact(score, 0.33401529220865844, diff=1e-3)


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
