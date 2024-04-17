import numpy as np

import scib
from tests.common import LOGGER, assert_near_exact


def _random_embedding(partition):
    """
    Copied from openproblems
    https://github.com/scottgigante-immunai/openproblems/blob/1fb510ed356925251cb0b9e11fbc7dcea5201161/openproblems/tasks/_batch_integration/batch_integration_embed/methods/baseline.py
    :param partition:
    :return:
    """
    from sklearn.preprocessing import LabelEncoder, OneHotEncoder

    embedding = OneHotEncoder().fit_transform(
        LabelEncoder().fit_transform(partition)[:, None]
    )
    embedding = embedding + np.random.uniform(-0.001, 0.001, embedding.shape)
    # convert to numpy array
    embedding = np.asarray(embedding)
    return embedding


def test_isolated_labels_f1(adata_neighbors):
    score = scib.me.isolated_labels_f1(
        adata_neighbors,
        label_key="celltype",
        batch_key="batch",
        embed="X_pca",
        verbose=True,
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.5581395348837209, diff=1e-12)


def test_isolated_labels_f1_precomputed(adata_clustered):
    score = scib.me.isolated_labels_f1(
        adata_clustered,
        label_key="celltype",
        batch_key="batch",
        embed="X_pca",
        cluster_key="cluster",
        verbose=True,
    )
    assert "iso_label" not in adata_clustered.obs.columns

    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.5581395348837209, diff=1e-12)


def test_isolated_labels_asw(adata_pca):
    score = scib.me.isolated_labels_asw(
        adata_pca,
        label_key="celltype",
        batch_key="batch",
        embed="X_pca",
        verbose=True,
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.5703811198472977, diff=1e-3)


def test_isolated_labels_f1_perfect(adata_pca):
    adata_pca.obsm["X_emb"] = _random_embedding(partition=adata_pca.obs["celltype"])
    score = scib.me.isolated_labels_f1(
        adata_pca,
        label_key="celltype",
        batch_key="batch",
        embed="X_emb",
        verbose=True,
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 1, diff=1e-12)


def test_isolated_labels_asw_perfect(adata_pca):
    adata_pca.obsm["X_emb"] = _random_embedding(partition=adata_pca.obs["celltype"])
    score = scib.me.isolated_labels_asw(
        adata_pca,
        label_key="celltype",
        batch_key="batch",
        embed="X_emb",
        verbose=True,
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 1, diff=1e-2)
