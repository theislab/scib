import scib
from tests.common import assert_near_exact


def test_scanvi(adata_paul15_template):
    adata = scib.ig.scanvi(
        adata_paul15_template, batch="batch", labels="celltype", max_epochs=20
    )

    scib.pp.reduce_data(
        adata, n_top_genes=200, neighbors=True, use_rep="X_emb", pca=True, umap=False
    )

    score = scib.me.graph_connectivity(adata, label_key="celltype")
    assert_near_exact(score, 0.9834078129657216, 1e-2)
