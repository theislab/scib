import scib
from tests.common import LOGGER, assert_near_exact


def test_clisi_full(adata):
    score = scib.me.clisi_graph(
        adata,
        batch_key="batch",
        label_key="celltype",
        scale=True,
        type_="full",
        verbose=True,
    )

    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.975, diff=1e-2)


def test_clisi_embed(adata_neighbors):
    adata_neighbors.obsm["X_emb"] = adata_neighbors.obsm["X_pca"]
    score = scib.me.clisi_graph(
        adata_neighbors,
        batch_key="batch",
        label_key="celltype",
        scale=True,
        type_="embed",
        verbose=True,
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.982, diff=1e-2)


def test_clisi_knn(adata_neighbors):
    score = scib.me.clisi_graph(
        adata_neighbors,
        batch_key="batch",
        label_key="celltype",
        scale=True,
        type_="graph",
        verbose=True,
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.982, diff=1e-2)


def test_clisi_full_parallel(adata):
    # test parallel
    score = scib.me.clisi_graph(
        adata,
        batch_key="batch",
        label_key="celltype",
        scale=True,
        type_="full",
        n_cores=2,
        verbose=True,
    )

    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.975, diff=1e-1)


def test_clisi_embed_parallel(adata_neighbors):
    adata_neighbors.obsm["X_emb"] = adata_neighbors.obsm["X_pca"]
    score = scib.me.clisi_graph(
        adata_neighbors,
        batch_key="batch",
        label_key="celltype",
        scale=True,
        type_="embed",
        n_cores=2,
        verbose=True,
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.982, diff=1e-2)


def test_clisi_knn_parallel(adata_neighbors):
    score = scib.me.clisi_graph(
        adata_neighbors,
        batch_key="batch",
        label_key="celltype",
        scale=True,
        type_="graph",
        n_cores=2,
        verbose=True,
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.982, diff=1e-2)
