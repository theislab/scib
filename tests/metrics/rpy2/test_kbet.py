import scib
from tests.common import LOGGER, assert_near_exact


def test_kbet(adata_pca):
    score = scib.me.kBET(
        adata_pca, batch_key="batch", label_key="celltype", embed="X_pca"
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.556108994805538, diff=1e-02)
