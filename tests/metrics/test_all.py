import scanpy as sc

import scib
from tests.common import LOGGER


def test_fast(adata_neighbors):
    metrics_df = scib.me.metrics_fast(
        adata_neighbors,
        adata_neighbors,
        batch_key="batch",
        label_key="celltype",
        embed="X_pca",
    )

    for score in metrics_df:
        LOGGER.info(f"score: {score}")
        assert 0 <= score <= 1


def test_slim(adata_paul15):
    sc.pp.pca(adata_paul15)
    sc.pp.neighbors(adata_paul15)
    sc.tl.dpt(adata_paul15)

    metrics_df = scib.me.metrics_slim(
        adata_paul15,
        adata_paul15,
        batch_key="batch",
        label_key="celltype",
        embed="X_pca",
    )

    for score in metrics_df:
        LOGGER.info(f"score: {score}")
        assert 0 <= score <= 1


# def test_all(adata_paul15):
#    sc.pp.pca(adata_paul15)
#    sc.pp.neighbors(adata_paul15)
#    sc.tl.dpt(adata_paul15)
#
#    metrics_df = scib.me.metrics_all(
#        adata_paul15,
#        adata_paul15,
#        batch_key='batch',
#        label_key='celltype',
#        embed='X_pca',
#    )
#
#    for score in metrics_df:
#        LOGGER.info(f"score: {score}")
#        assert 0 <= score <= 1
