import scib
from tests.common import LOGGER, assert_near_exact


def test_isolated_labels_F1(adata_neighbors):
    score = scib.me.isolated_labels_f1(
        adata_neighbors,
        label_key="celltype",
        batch_key="batch",
        embed="X_pca",
        verbose=True,
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.5581395348837209, diff=1e-12)


def test_isolated_labels_ASW(adata_neighbors):
    score = scib.me.isolated_labels_asw(
        adata_neighbors,
        label_key="celltype",
        batch_key="batch",
        embed="X_pca",
        verbose=True,
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.6101431176066399, diff=1e-3)
