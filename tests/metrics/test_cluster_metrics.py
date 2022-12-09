import scib
from tests.common import LOGGER, assert_near_exact


def test_nmi_trivial(adata):
    score = scib.me.nmi(adata, cluster_key="celltype", label_key="celltype")
    assert_near_exact(score, 1, diff=1e-12)


def test_nmi_deprecated(adata):
    scib.me.nmi(adata, cluster_key="celltype", label_key="celltype")


def test_ari_trivial(adata):
    score = scib.me.ari(adata, cluster_key="celltype", label_key="celltype")
    assert_near_exact(score, 1, diff=1e-12)


def test_ari_deprecated(adata):
    scib.me.ari(adata, group1="celltype", group2="celltype")


def test_nmi(adata_neighbors):
    _, _, nmi_all = scib.cl.opt_louvain(
        adata_neighbors,
        label_key="celltype",
        cluster_key="cluster",
        function=scib.me.nmi,
        plot=False,
        inplace=True,
        force=True,
        verbose=True,
    )

    for score in nmi_all["score"]:
        assert 0 <= score <= 1


def test_ari(adata_clustered):
    score = scib.me.ari(adata_clustered, cluster_key="cluster", label_key="celltype")
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.7614422905830917, diff=1e-2)
