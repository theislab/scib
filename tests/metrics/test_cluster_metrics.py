import scib
from tests.common import LOGGER, assert_near_exact


def test_deprecated(adata):
    scib.me.nmi(adata, group1="celltype", group2="celltype")
    scib.me.ari(adata, group1="celltype", group2="celltype")


def test_nmi_trivial(adata):
    score = scib.me.nmi(adata, cluster_key="celltype", label_key="celltype")
    assert_near_exact(score, 1, diff=1e-12)


def test_ari_trivial(adata):
    score = scib.me.ari(adata, cluster_key="celltype", label_key="celltype")
    assert_near_exact(score, 1, diff=1e-12)


def test_nmi(adata_clustered):
    score = scib.me.nmi(adata_clustered, cluster_key="cluster", label_key="celltype")
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.7424614219722735, diff=1e-2)


def test_ari(adata_clustered):
    score = scib.me.ari(adata_clustered, cluster_key="cluster", label_key="celltype")
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.7614422905830917, diff=1e-2)
