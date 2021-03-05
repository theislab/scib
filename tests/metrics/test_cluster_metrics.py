from tests.common import *


def test_nmi_trivial(adata):
    score = scIB.me.nmi(adata, 'celltype', 'celltype')
    assert score == 1


def test_ari_trivial(adata):
    score = scIB.me.ari(adata, 'celltype', 'celltype')
    assert score == 1


def test_nmi(adata, adata_neighbors, cluster_factory):
    _, _, nmi_all = cluster_factory(
        adata_neighbors(adata),
        cluster_key='cluster',
        label_key='celltype',
        verbose=True
    )
    for score in nmi_all['score']:
        assert score >= 0
        assert score <= 1


def test_ari(adata, adata_clustered):
    score = scIB.me.ari(adata_clustered(adata), 'cluster', 'celltype')
    LOGGER.info(f"score: {score}")
    assert score >= 0
    assert score <= 1


def test_isolated_labels_F1(adata, adata_neighbors):
    score = scIB.me.isolated_labels(
        adata_neighbors(adata),
        label_key='celltype',
        batch_key='batch',
        embed='X_pca',
        cluster=True,
        verbose=True
    )
    LOGGER.info(f"score: {score}")
    assert score <= 1
    assert score >= 0
