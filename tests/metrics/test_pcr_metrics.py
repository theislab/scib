from tests.common import *


def test_pc_regression(adata):
    scib.me.pcr.pc_regression(adata.X, adata.obs["batch"])


def test_pcr_batch(adata):
    # no PCA precomputed
    score = scib.me.pcr_comparison(
        adata, adata,
        covariate='batch',
        n_comps=50,
        scale=True
    )
    LOGGER.info(f"no PCA precomputed: {score}")
    assert_near_exact(score, 0, diff=1e-6)


def test_pcr_batch_precomputed(adata_pca):
    score = scib.me.pcr_comparison(adata_pca, adata_pca, covariate='batch', scale=True)
    LOGGER.info(f"precomputed PCA: {score}")
    assert_near_exact(score, 0, diff=1e-6)


def test_pcr_batch_embedding(adata):
    # use different embedding
    score = scib.me.pcr_comparison(
        adata_pre=adata,
        adata_post=add_embed(adata, type_='full'),
        covariate='batch',
        embed='X_emb',
        n_comps=50,
        scale=True
    )
    LOGGER.info(f"using embedding: {score}")
    assert_near_exact(score, 0, diff=1e-6)
