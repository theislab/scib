import scib
from tests.common import assert_near_exact


def test_scvi(adata_paul15_template):
    adata = scib.ig.scvi(adata_paul15_template, batch="batch", max_epochs=20)

    scib.pp.reduce_data(
        adata, n_top_genes=200, neighbors=True, use_rep="X_emb", pca=True, umap=False
    )

    score = scib.me.graph_connectivity(adata, label_key="celltype")
    assert_near_exact(score, 0.9684638088694193, 1e-2)
