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


def test_metrics_uses_embed_for_clisi():
    """Regression test: metrics() must pass embed through to LISI scores.

    Build a small AnnData that has the embed rep under a non-default obsm
    key and no X_emb. Before the fix, clisi_graph/ilisi_graph fell back to
    the X_emb default and raised a KeyError.
    """
    import anndata as ad
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(0)
    n = 60
    adata = ad.AnnData(
        X=rng.normal(size=(n, 10)),
        obs=pd.DataFrame(
            {
                "batch": ["a", "b"] * (n // 2),
                "celltype": ["x", "y"] * (n // 2),
            },
            index=[f"cell{i}" for i in range(n)],
        ),
    )
    adata.obs["batch"] = adata.obs["batch"].astype("category")
    adata.obs["celltype"] = adata.obs["celltype"].astype("category")
    adata.obsm["X_custom"] = rng.normal(size=(n, 5))
    assert "X_emb" not in adata.obsm

    scores = scib.me.metrics(
        adata,
        adata,
        batch_key="batch",
        label_key="celltype",
        embed="X_custom",
        type_="embed",
        clisi_=True,
        ilisi_=True,
    )
    assert np.isfinite(scores.loc["cLISI", 0])
    assert np.isfinite(scores.loc["iLISI", 0])
