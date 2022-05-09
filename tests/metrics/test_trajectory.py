import scanpy as sc

import scib
from tests.common import LOGGER, assert_near_exact


def test_trajectory(adata_neighbors):
    adata_int = adata_neighbors.copy()

    adata_neighbors.uns["iroot"] = 0
    sc.tl.diffmap(adata_neighbors)
    sc.tl.dpt(adata_neighbors)

    score = scib.me.trajectory_conservation(
        adata_pre=adata_neighbors,
        adata_post=adata_int,
        label_key="celltype",
        pseudotime_key="dpt_pseudotime",
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.9561, diff=1e-5)


def test_trajectory_batch(adata_neighbors):
    adata_int = adata_neighbors.copy()

    adata_neighbors.uns["iroot"] = 0
    sc.tl.diffmap(adata_neighbors)
    sc.tl.dpt(adata_neighbors)

    score = scib.me.trajectory_conservation(
        adata_pre=adata_neighbors,
        adata_post=adata_int,
        label_key="celltype",
        batch_key="batch",
        pseudotime_key="dpt_pseudotime",
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.96317, diff=1e-5)
