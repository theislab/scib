import pandas as pd

import scib
from tests.common import LOGGER, assert_near_exact


def test_scanorama(adata_paul15_template):
    adata = scib.ig.scanorama(adata_paul15_template, batch="batch")

    # check full feature output
    scib.pp.reduce_data(
        adata, pca=True, n_top_genes=200, neighbors=True, use_rep="X_pca", umap=False
    )
    score = scib.me.graph_connectivity(adata, label_key="celltype")
    LOGGER.info(f"\nscore: {score}")

    # check embedding output
    scib.pp.reduce_data(adata, pca=False, neighbors=True, use_rep="X_emb", umap=False)
    score = scib.me.graph_connectivity(adata, label_key="celltype")
    LOGGER.info(f"\nscore: {score}")

    # check NMI after clustering
    res_max, score_max, _ = scib.cl.opt_louvain(
        adata,
        label_key="celltype",
        cluster_key="cluster",
        plot=False,
        inplace=True,
    )
    LOGGER.info(f"max resolution: {res_max}, max score: {score_max}")

    assert_near_exact(score_max, 0.6440883682371078, 1e-2)


def test_scanorama_batch_cols(adata_paul15_template):
    adata = scib.ig.scanorama(adata_paul15_template, batch="batch")

    # check that batches match before and after integration
    batch_before = adata_paul15_template.obs["batch"].value_counts()
    batch_after = adata.obs["batch"].value_counts()

    pd.testing.assert_series_equal(batch_before, batch_after)


# def test_scanorama_dup_cols(adata_paul15_template):
#     adata = adata_paul15_template
#     adata.var[['1', '2']] = (0, 2)
#     # create adata with multiple variables with same name
#     adata.var = adata.var.rename(
#         columns={name: 'ambig_column' for name in adata.var.columns.values}
#     )
#
#     LOGGER.info(adata)
#
#     # concatenate adatas with ambiguous column naming
#     anndata.AnnData.concatenate(
#         adata, adata, index_unique=None
#     )
