import scib
from tests.common import LOGGER, assert_near_exact


def test_scanvi(adata_paul15_template):
    adata = scib.ig.scanvi(
        adata_paul15_template, batch="batch", labels="celltype", max_epochs=20
    )

    scib.pp.reduce_data(
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
    LOGGER.info(f"max resolution: {res_max}, max score: {score_max}")

    assert_near_exact(score_max, 0.8679758761826293, 1e-2)
