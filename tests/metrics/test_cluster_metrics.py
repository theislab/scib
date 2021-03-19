from tests.common import *


def test_nmi_trivial(adata):
    score = scIB.me.nmi(adata, 'celltype', 'celltype')
    assert score == 1


def test_ari_trivial(adata):
    score = scIB.me.ari(adata, 'celltype', 'celltype')
    assert score == 1


def test_nmi(adata_neighbors):

    _, _, nmi_all = scIB.cl.opt_louvain(
        adata_neighbors,
        label_key='celltype',
        cluster_key='cluster',
        function=scIB.me.nmi,
        plot=False,
        inplace=True,
        force=True,
        verbose=True
    )

    for score in nmi_all['score']:
        assert 0 <= score <= 1


def test_ari(adata_clustered):
    score = scIB.me.ari(adata_clustered, group1='cluster', group2='celltype')
    LOGGER.info(f"score: {score}")
    assert 0 <= score <= 1


def test_isolated_labels_F1(adata_neighbors):
    score = scIB.me.isolated_labels(
        adata_neighbors,
        label_key='celltype',
        batch_key='batch',
        embed='X_pca',
        cluster=True,
        verbose=True
    )
    LOGGER.info(f"score: {score}")
    assert 0 <= score <= 1
