from tests.common import *


def test_lisi_full(adata):
    ilisi, clisi = scIB.me.lisi_graph(
        adata,
        batch_key='batch',
        label_key='celltype',
        scale=True,
        type_='full'
    )

    LOGGER.info(f"score: {ilisi, clisi}")
    assert 0 <= ilisi <= 1
    assert 0 <= clisi <= 1


def test_lisi_embed(adata_neighbors):
    adata_neighbors.obsm['X_emb'] = adata_neighbors.obsm['X_pca']
    ilisi, clisi = scIB.me.lisi_graph(
        adata_neighbors,
        batch_key='batch',
        label_key='celltype',
        scale=True,
        type_='embed'
    )
    LOGGER.info(f"score: {ilisi, clisi}")
    assert 0 <= ilisi <= 1
    assert 0 <= clisi <= 1


def test_lisi_graph(adata_neighbors):
    ilisi, clisi = scIB.me.lisi_graph(
        adata_neighbors,
        batch_key='batch',
        label_key='celltype',
        scale=True,
        type_='graph',

    )
    LOGGER.info(f"score: {ilisi, clisi}")
    assert 0 <= ilisi <= 1
    assert 0 <= clisi <= 1
