import pandas as pd

import scib
from tests.common import LOGGER, assert_near_exact


def test_scanorama(adata_paul15_template):
    adata = scib.ig.scanorama(adata_paul15_template, batch="batch")

    scib.pp.reduce_data(
        adata, n_top_genes=200, neighbors=True, use_rep="X_emb", pca=True, umap=False
    )

    # check NMI after clustering
    res_max, score_max, _ = scib.cl.cluster_optimal_resolution(
        adata, label_key="celltype", cluster_key="cluster", return_all=True
    )
    LOGGER.info(f"max resolution: {res_max}, max score: {score_max}")

    assert_near_exact(score_max, 0.6610082444492823, 1e-2)


def test_scanorama_batch_cols(adata_paul15_template):
    adata = scib.ig.scanorama(adata_paul15_template, batch="batch")

    # check that batches match before and after integration
    batch_before = adata_paul15_template.obs["batch"].value_counts()
    batch_after = adata.obs["batch"].value_counts()

    pd.testing.assert_series_equal(batch_before, batch_after)
