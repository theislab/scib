from tests.common import *


def test_cell_cycle(adata_paul15):
    adata = adata_paul15
    adata_int = adata.copy()

    # only final score
    score = scIB.me.cell_cycle(
        adata,
        adata_int,
        batch_key='batch',
        organism='mouse',
        #recompute_cc=True,
        verbose=True
    )
    LOGGER.info(f"score: {score}")
    assert score == 1


def test_cell_cycle_all(adata_paul15):
    adata = adata_paul15
    adata_int = adata.copy()

    # get all intermediate scores
    scores_df = scIB.me.cell_cycle(
        adata,
        adata_int,
        batch_key='batch',
        organism='mouse',
        #recompute_cc=True,
        agg_func=None,
        verbose=True
    )
    LOGGER.info(f"\nscore: {scores_df}")
    assert isinstance(scores_df, pd.DataFrame)
    for i in scores_df['score']:
        assert i == 1


def test_hvg_overlap(adata):
    adata_int = adata.copy()
    score = scIB.me.hvg_overlap(
        adata_int,
        adata,
        batch='batch',
        n_hvg=500
    )
    LOGGER.info(f"score: {score}")
    assert score == 1


# TODO: trajectory conservation metric
