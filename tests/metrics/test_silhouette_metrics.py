import numpy as np

import scib
from tests.common import LOGGER, assert_near_exact


def test_silhouette(adata_pca):
    score = scib.me.silhouette(
        adata_pca, label_key="celltype", embed="X_pca", scale=True
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.5626532882452011, diff=1e-2)


def test_silhouette_batch(adata_pca):
    score = scib.me.silhouette_batch(
        adata_pca,
        batch_key="batch",
        label_key="celltype",
        embed="X_pca",
        scale=True,
        verbose=False,
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.9014384369842835, diff=1e-2)


def test_silhouette_batch_empty(adata_pca):
    adata_pca.obs["label"] = "label"
    adata_pca.obs["batch"] = "batch"
    asw, sil_means, sil_df = scib.me.silhouette_batch(
        adata_pca,
        batch_key="batch",
        label_key="celltype",
        embed="X_pca",
        scale=True,
        verbose=False,
        return_all=True,
    )
    assert np.isnan(asw)
    assert np.isnan(sil_means)
    assert sil_df.shape[0] == 0
    LOGGER.info(sil_df)
