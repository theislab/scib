from tests.common import *


def test_trajectory(adata_neighbors):
    sc.tl.dpt(adata_neighbors)
    adata_int = adata_neighbors.copy()
    score = scib.me.hvg_overlap(
        adata_int,
        adata_neighbors,
        batch='batch',
        n_hvg=500
    )
    LOGGER.info(f"score: {score}")
    assert score == 1
