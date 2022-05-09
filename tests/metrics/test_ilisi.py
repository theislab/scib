import scib
from tests.common import LOGGER, assert_near_exact


def test_ilisi_full(adata):
    score = scib.me.ilisi_graph(
        adata, batch_key="batch", scale=True, type_="full", verbose=True
    )

    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.235, diff=1e-2)


def test_ilisi_embed(adata_neighbors):
    adata_neighbors.obsm["X_emb"] = adata_neighbors.obsm["X_pca"]
    score = scib.me.ilisi_graph(
        adata_neighbors, batch_key="batch", scale=True, type_="embed", verbose=True
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.238, diff=1e-2)


def test_ilisi_knn(adata_neighbors):
    score = scib.me.ilisi_graph(
        adata_neighbors, batch_key="batch", scale=True, type_="graph", verbose=True
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.238, diff=1e-2)


def test_ilisi_full_parallel(adata):
    score = scib.me.ilisi_graph(
        adata, batch_key="batch", scale=True, type_="full", verbose=True
    )

    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.235, diff=1e-2)


def test_ilisi_embed_parallel(adata_neighbors):
    adata_neighbors.obsm["X_emb"] = adata_neighbors.obsm["X_pca"]
    score = scib.me.ilisi_graph(
        adata_neighbors, batch_key="batch", scale=True, type_="embed", verbose=True
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.238, diff=1e-2)


def test_ilisi_knn_parallel(adata_neighbors):
    score = scib.me.ilisi_graph(
        adata_neighbors, batch_key="batch", scale=True, type_="graph", verbose=True
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.238, diff=1e-2)
