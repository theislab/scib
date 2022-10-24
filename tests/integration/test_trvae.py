import scib
from tests.common import LOGGER


def test_trvae(adata_paul15_template):
    adata = scib.ig.trvae(adata_paul15_template, batch="batch")

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
    # assert_near_exact(score, 0.9834078129657216, 1e-2)
