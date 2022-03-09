from tests.common import *


def test_nmi_trivial(adata):
    score = scib.me.nmi(adata, 'celltype', 'celltype')
    assert_near_exact(score, 1, diff=1e-12)


def test_ari_trivial(adata):
    score = scib.me.ari(adata, 'celltype', 'celltype')
    assert_near_exact(score, 1, diff=1e-12)


def test_nmi(adata_neighbors):
    _, _, nmi_all = scib.cl.opt_louvain(
        adata_neighbors,
        label_key='celltype',
        cluster_key='cluster',
        function=scib.me.nmi,
        plot=False,
        inplace=True,
        force=True,
        verbose=True
    )

    for score in nmi_all['score']:
        assert 0 <= score <= 1


def test_ari(adata_clustered):
    score = scib.me.ari(adata_clustered, group1='cluster', group2='celltype')
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.7614422905830917, diff=1e-2)


def test_isolated_labels_F1(adata_neighbors):
    score = scib.me.isolated_labels(
        adata_neighbors,
        label_key='celltype',
        batch_key='batch',
        embed='X_pca',
        cluster=True,
        verbose=True
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.5581395348837209, diff=1e-12)
