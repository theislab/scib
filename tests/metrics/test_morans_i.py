import scib
from tests.common import LOGGER, assert_near_exact


def test_morans_i_nhvg(adata_neighbors):
    adata_int = adata_neighbors.copy()
    score = scib.me.morans_i(
        adata_pre=adata_neighbors,
        adata_post=adata_int,
        batch_key="batch",
        n_hvg=100,
    )
    LOGGER.info(f"score: {score}")
    expected_score = 0.65204
    assert_near_exact(score, expected_score, diff=1e-5)


def test_morans_i_hvgs(adata_neighbors):
    adata_int = adata_neighbors.copy()
    score = scib.me.morans_i(
        adata_pre=adata_neighbors,
        adata_post=adata_int,
        batch_key="batch",
        hvgs=[
            "TNFRSF4",
            "TNFRSF1B",
            "EFHD2",
            "C1QA",
            "C1QB",
            "STMN1",
            "SMAP2",
            "PRDX1",
            "SCP2",
            "C1orf162",
        ],
    )
    LOGGER.info(f"score: {score}")
    expected_score = 0.64337
    assert_near_exact(score, expected_score, diff=1e-5)
