import pandas as pd
from scipy.sparse import csr_matrix

import scib
from tests.common import LOGGER, assert_near_exact


def test_cell_cycle(adata_paul15):
    adata = adata_paul15
    adata_int = adata.copy()

    # only final score
    score = scib.me.cell_cycle(
        adata,
        adata_int,
        batch_key="batch",
        organism="mouse",
        # recompute_cc=True,
        verbose=True,
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 1, diff=1e-12)


def test_cell_cycle_sparse(adata_paul15):
    adata = adata_paul15
    adata_int = adata.copy()

    # sparse matrix
    adata.X = csr_matrix(adata.X)
    adata_int.X = csr_matrix(adata.X)

    # only final score
    score = scib.me.cell_cycle(
        adata,
        adata_int,
        batch_key="batch",
        organism="mouse",
        n_comps=adata.shape[1],
        verbose=True,
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 1, diff=1e-12)


def test_cell_cycle_all(adata_paul15):
    adata = adata_paul15
    adata_int = adata.copy()

    # get all intermediate scores
    scores_df = scib.me.cell_cycle(
        adata,
        adata_int,
        batch_key="batch",
        organism="mouse",
        # recompute_cc=True,
        agg_func=None,
        verbose=True,
    )
    LOGGER.info(f"\nscore: {scores_df}")
    assert isinstance(scores_df, pd.DataFrame)
    for i in scores_df["score"]:
        assert_near_exact(i, 1, diff=1e-12)


def test_hvg_overlap(adata):
    adata_int = adata.copy()
    score = scib.me.hvg_overlap(adata_int, adata, batch_key="batch", n_hvg=500)
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 1, diff=1e-12)
