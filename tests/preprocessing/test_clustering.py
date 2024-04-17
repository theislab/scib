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
        resolutions=scib.cl.get_resolutions(n=10, min=0, max=1),
        return_all=True,
    )
    assert isinstance(score_all, pd.DataFrame)
    assert "cluster" in adata_neighbors.obs.columns

    LOGGER.info(f"max resolution: {res_max}, max score: {score_max}")
    assert res_max == 0.7
    assert_near_exact(score_max, 0.7432787576640969, diff=1e-3)


def test_cluster_optimal_resolution_leiden(adata_neighbors):
    res_max, score_max, score_all = scib.cl.cluster_optimal_resolution(
        adata_neighbors,
        label_key="celltype",
        cluster_key="cluster",
        cluster_function=sc.tl.leiden,
        resolutions=scib.cl.get_resolutions(n=10, min=0, max=1),
        return_all=True,
    )
    assert isinstance(score_all, pd.DataFrame)
    assert "cluster" in adata_neighbors.obs.columns

    LOGGER.info(f"max resolution: {res_max}, max score: {score_max}")
    assert res_max == 0.5
    assert_near_exact(score_max, 0.7424614219722735, diff=1e-3)


def test_precomputed_cluster(adata):
    resolutions = scib.cl.get_resolutions(n=10, min=0, max=1)
    for res in resolutions:
        adata.obs[f"cluster_{res}"] = adata.obs["celltype"]

    res_max, score_max = scib.cl.cluster_optimal_resolution(
        adata,
        cluster_key="cluster",
        label_key="celltype",
        force=False,
        resolutions=resolutions,
    )
    assert res_max == 0.1
    assert_near_exact(score_max, 1, diff=0)


def test_precomputed_cluster_force(adata_neighbors):
    resolutions = scib.cl.get_resolutions(n=10, min=0, max=1)
    for res in resolutions:
        adata_neighbors.obs[f"cluster_{res}"] = adata_neighbors.obs["celltype"]

    res_max, score_max = scib.cl.cluster_optimal_resolution(
        adata_neighbors,
        cluster_key="cluster",
        label_key="celltype",
        resolutions=resolutions,
        force=True,
    )
    assert res_max == 0.5
    assert_near_exact(score_max, 0.7424614219722736, diff=1e-5)
