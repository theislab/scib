import scib
from tests.common import LOGGER, assert_near_exact


def test_harmony(adata_paul15_template):
    adata = scib.ig.harmony(adata_paul15_template, batch="batch")

    scib.pp.reduce_data(
        adata, n_top_genes=200, neighbors=True, use_rep="X_emb", pca=True, umap=False
    )

    # check NMI after clustering
    res_max, score_max, _ = scib.cl.cluster_optimal_resolution(
        adata,
        label_key="celltype",
        cluster_key="cluster",
        use_rep="X_emb",
        return_all=True,
    )
    LOGGER.info(f"max resolution: {res_max}, max score: {score_max}")

    assert_near_exact(score_max, 0.5050152860428626, 1e-2)
