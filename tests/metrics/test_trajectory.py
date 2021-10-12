from tests.common import *


def test_trajectory(adata_neighbors):
    adata_int = adata_neighbors.copy()

    adata_neighbors.uns['iroot'] = 0
    sc.tl.diffmap(adata_neighbors)
    sc.tl.dpt(adata_neighbors)

    score = scIB.me.trajectory_conservation(
        adata_pre=adata_neighbors,
        adata_post=adata_int,
        label_key='celltype',
        pseudotime_key='dpt_pseudotime'
    )
    LOGGER.info(f"score: {score}")
    assert 0.956091 <= score <= 0.956092


def test_trajectory_batch(adata_neighbors):
    adata_int = adata_neighbors.copy()

    adata_neighbors.uns['iroot'] = 0
    sc.tl.diffmap(adata_neighbors)
    sc.tl.dpt(adata_neighbors)

    score = scIB.me.trajectory_conservation(
        adata_pre=adata_neighbors,
        adata_post=adata_int,
        label_key='celltype',
        batch_key='batch',
        pseudotime_key='dpt_pseudotime'
    )
    LOGGER.info(f"score: {score}")
    assert 0.9631630 <= score <= 0.9631631
