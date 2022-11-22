import pandas as pd
import scanpy as sc

import scib
from tests.common import LOGGER, assert_near_exact


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
    assert_near_exact(score_max, 0.7432787576640969, diff=1e-3)


def test_cluster_optimal_resolution_louvain(adata_neighbors):
    res_max, score_max, score_all = scib.cl.cluster_optimal_resolution(
        adata_neighbors,
        label_key="celltype",
        cluster_key="cluster",
        cluster_function=sc.tl.louvain,
        resolutions=scib.cl.get_resolutions(n=20, min=0.1, max=2),
        return_all=True,
    )
    assert isinstance(score_all, pd.DataFrame)

    LOGGER.info(f"max resolution: {res_max}, max score: {score_max}")
    assert res_max == 0.7
    assert_near_exact(score_max, 0.7432787576640969, diff=1e-3)


def test_cluster_optimal_resolution_leiden(adata_neighbors):
    res_max, score_max, score_all = scib.cl.cluster_optimal_resolution(
        adata_neighbors,
        label_key="celltype",
        cluster_key="cluster",
        cluster_function=sc.tl.leiden,
        return_all=True,
    )
    assert isinstance(score_all, pd.DataFrame)

    LOGGER.info(f"max resolution: {res_max}, max score: {score_max}")
    assert res_max == 0.5
    assert_near_exact(score_max, 0.7424614219722735, diff=1e-3)
