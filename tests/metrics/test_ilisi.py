from tests.common import *


def test_ilisi_full(adata):
    score = scIB.me.ilisi_graph(
        adata,
        batch_key='batch',
        scale=True,
        type_='full'
    )

    LOGGER.info(f"score: {score}")
    assert 0.234 <= score <= 0.2341


def test_ilisi_embed(adata_neighbors):
    adata_neighbors.obsm['X_emb'] = adata_neighbors.obsm['X_pca']
    score = scIB.me.ilisi_graph(
        adata_neighbors,
        batch_key='batch',
        scale=True,
        type_='embed'
    )
    LOGGER.info(f"score: {score}")
    assert 0.2377 <= score <= 0.2378


def test_ilisi_knn(adata_neighbors):
    score = scIB.me.ilisi_graph(
        adata_neighbors,
        batch_key='batch',
        scale=True,
        type_='graph'
    )
    LOGGER.info(f"score: {score}")
    assert 0.2377 <= score <= 0.2378
