from tests.common import *


def test_ilisi_full(adata):
    score = scib.me.ilisi_graph(
        adata,
        batch_key='batch',
        scale=True,
        type_='full'
    )

    LOGGER.info(f"score: {score}")
    assert 0.234 <= score <= 0.235


def test_ilisi_embed(adata_neighbors):
    adata_neighbors.obsm['X_emb'] = adata_neighbors.obsm['X_pca']
    score = scib.me.ilisi_graph(
        adata_neighbors,
        batch_key='batch',
        scale=True,
        type_='embed'
    )
    LOGGER.info(f"score: {score}")
    assert 0.237 <= score <= 0.238


def test_ilisi_knn(adata_neighbors):
    score = scib.me.ilisi_graph(
        adata_neighbors,
        batch_key='batch',
        scale=True,
        type_='graph'
    )
    LOGGER.info(f"score: {score}")
    assert 0.237 <= score <= 0.238
