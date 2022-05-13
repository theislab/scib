import pandas as pd

import scib
from tests.common import LOGGER, assert_near_exact


def test_scanorama(adata_paul15_template):
    adata = scib.ig.scanorama(adata_paul15_template, batch="batch")

    # check that batches match before and after integration
    batch_before = adata_paul15_template.obs["batch"].value_counts()
    batch_after = adata.obs["batch"].value_counts()

    pd.testing.assert_series_equal(batch_before, batch_after)

    scib.preprocessing.reduce_data(
        adata, n_top_genes=200, neighbors=True, use_rep="X_emb", pca=True, umap=False
    )

    # check NMI after clustering
    res_max, score_max, _ = scib.cl.opt_louvain(
        adata,
        label_key="celltype",
        cluster_key="cluster",
        plot=False,
        inplace=True,
    )
    LOGGER.info(f"max resolution: {res_max}, max NMI: {score_max}")

    assert res_max == 2.0
    assert_near_exact(score_max, 0.6440883682371078, 1e-2)
