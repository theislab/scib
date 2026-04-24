import scib
from tests.common import LOGGER, assert_near_exact


def test_scvi(adata_paul15_template):
    adata = scib.ig.scvi(adata_paul15_template, batch="batch", max_epochs=20)

    scib.pp.reduce_data(
        adata, n_top_genes=200, neighbors=True, use_rep="X_emb", pca=True, umap=False
    )

    score = scib.me.graph_connectivity(adata, label_key="celltype")
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.96, 1e-1)
