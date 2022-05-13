import pandas as pd

import scib
from tests.common import LOGGER


def test_opt_louvain(adata_neighbors):
    res_max, score_max, score_all, clustering = scib.cl.opt_louvain(
        adata_neighbors,
        label_key="celltype",
        cluster_key="cluster",
        plot=False,
        inplace=False,
    )
    assert isinstance(score_all, pd.DataFrame)
    assert isinstance(clustering, pd.Series)

    LOGGER.info(f"max resolution: {res_max}, max score: {score_max}")
    assert res_max == 0.7
    assert score_max == 0.7432787576640969
