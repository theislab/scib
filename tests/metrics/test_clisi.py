from tests.common import *


def test_clisi_full(adata):
    score = scib.me.clisi_graph(
        adata,
        batch_key='batch',
        label_key='celltype',
        scale=True,
        type_='full'
    )

    LOGGER.info(f"score: {score}")
    assert 0.822 <= score <= 0.823


def test_clisi_embed(adata_neighbors):
    adata_neighbors.obsm['X_emb'] = adata_neighbors.obsm['X_pca']
    score = scib.me.clisi_graph(
        adata_neighbors,
        batch_key='batch',
        label_key='celltype',
        scale=True,
        type_='embed'
    )
    LOGGER.info(f"score: {score}")
    assert 0.868 <= score <= 0.869


def test_clisi_knn(adata_neighbors):
    score = scib.me.clisi_graph(
        adata_neighbors,
        batch_key='batch',
        label_key='celltype',
        scale=True,
        type_='graph'
    )
    LOGGER.info(f"score: {score}")
    assert 0.868 <= score <= 0.869
